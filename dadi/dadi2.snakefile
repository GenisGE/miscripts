"""
  snakemake to run dadi. config will need:
    pairs: dic of population pairs to estimate demographic inference, key:values are poppair:[pop1, pop2]
    models: dic of dic of dics, second level are model names, then for each model a dic with
           funcname: name of function specifying model; 
           model_params: list of parameter names (MUST be pop sizes start with N, times with T, migration rates with m) (if we want to include other kind of params, like fractions, would need to change reescaling rule)
            upper_bound: list with parameter upper_bounds
            lower_bound: list with parameter lower bounds
    popn: dic with sample sizes per population, key value are pop:N where N is an integer
    outmain: path to main otuput oflder
    nruns: fixed number of independent optimization runs to do first, that are ran in parallel and then checked for convergende
    maxn: maximum number of optimization runs to do if not converging
    mu: per generation mutation rate
    g: generation time

 optimization rutine first does a fixed number of runs per model (nruns, e.g. 10). Then checks for convergence and if not converged keeps doing runs until convergence or reaches maximum number without convergence (maxn, e.g. 50). Convergence defined heuristically as having at least 3 top likelihood runs meeting (abs(delta_like) < 4).

This one will use the perturbe parameters to chose initital guess for each run, rahter than sampling over whole parameter space that might be causing very slow convergence.
"""

import dadi
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import sys
sys.path.append("scripts")
import dadi_models

OUTMAIN=config["outmain"]

# for now import 2dsfs from fst workflow. Ideally in future create sfs only snakemake from which both fst and dadi would import 2dsfs
subworkflow fst:
    workdir:
        "../fst/"
    snakefile:
        "../fst/fst.snakefile"
    configfile:
        "../fst/config2.yaml"

rule all:
    input:
        expand(os.path.join(OUTMAIN, "2dsfs_dadi", "{p[0]}_{p[1]}_{w}_dadi.png"), p = config["pairs"].values(), w=["folded", "unfolded"]),
        expand(os.path.join(OUTMAIN, "fit", "{pair}", "{model}", "maxlike_params.txt"),
               pair=config["pairs"].keys(), model = config["models"].keys()),
        expand(os.path.join(OUTMAIN, "fit", "{pair}", "{model}", "collected_pars_allruns.txt"),
               pair=config["pairs"].keys(), model = config["models"].keys()),
        expand(os.path.join(OUTMAIN, "fit", "{pair}", "{model}", "{pair}_{model}_dadifit_plot.png"),
               pair=config["pairs"].keys(), model = config["models"].keys()),
        expand(os.path.join(OUTMAIN, "fit", "{pair}", "{model}", "rescaled_maxlike_params.txt"),
               pair=config["pairs"].keys(), model = config["models"].keys()),
        expand(os.path.join(OUTMAIN, "fit", "{pair}", "model_selection_aic.tsv"),
               pair=config["pairs"].keys(), model= config["models"].keys())





rule dadi_format_sfs:
    input:
        sfs = fst(os.path.join("results2", "2dsfs", "{p1}_{p2}.sfs"))
    output:
        sfs_unfolded = os.path.join(OUTMAIN, "2dsfs_dadi", "{p1}_{p2}_unfolded_dadi.sfs"),
#        sfs_folded = os.path.join(OUTMAIN, "2dsfs_dadi", "{p1}_{p2}_folded_dadi.sfs")
    params:
        pop1 = "{p1}",
        pop2 = "{p2}",
        n1 = lambda wildcards: config["popn"][wildcards.p1],
        n2 = lambda wildcards: config["popn"][wildcards.p2]
    run:
        comment = f'#unfolded sfs {params.pop1} {params.pop2}\n'
        n1 = int(params.n1) * 2 + 1
        n2 = int(params.n2) * 2 + 1
        header = f'{n1} {n2} unfolded\n'
        with open(input.sfs, "r") as fh:
            sfs = fh.read()
        mask = "1 " + "0 " * (len(sfs)-2) + " 1\n" # SHOULD CHECK IF WE REALLY WANT TO MASK FIXED SITES WHEN USING SEQUENCING DATA
        with open(output.sfs_unfolded, "w+") as fh:
            fh.write(comment)
            fh.write(header)
            fh.write(sfs)
            fh.write(mask)



rule plot_dadi_sfs:
    input:
        sfs = rules.dadi_format_sfs.output.sfs_unfolded
    output:
        sfs_unfolded_plot = os.path.join(OUTMAIN, "2dsfs_dadi", "{p1}_{p2}_unfolded_dadi.png"),
        sfs_folded_plot = os.path.join(OUTMAIN, "2dsfs_dadi", "{p1}_{p2}_folded_dadi.png")
    params:
        pop1 = "{p1}",
        pop2 = "{p2}"
    run:
        fs = dadi.Spectrum.from_file(input.sfs)
        fs.pop_ids = [params.pop1, params.pop2]
        dadi.Plotting.plot_single_2d_sfs(fs, vmin=1)
        plt.savefig(output.sfs_unfolded_plot)
        plt.clf()
        fs_fold = fs.fold()
        dadi.Plotting.plot_single_2d_sfs(fs_fold, vmin=1)
        plt.savefig(output.sfs_folded_plot)
        plt.clf()




rule fit_model:
    input:
        dadisfs = lambda wildcards: os.path.join(OUTMAIN, "2dsfs_dadi", "{p[0]}_{p[1]}_unfolded_dadi.sfs".format(p=config["pairs"][wildcards.pair])), # NOT SURE THIS WILL WORK
    output:
        maxlpars = temp(os.path.join(OUTMAIN, "fit", "{pair}", "{model}", "run_{seed}.maxlpars")),
        like = temp(os.path.join(OUTMAIN, "fit", "{pair}", "{model}", "run_{seed}.likehood"))
    params:
        seed = "{seed}",
        funcname = lambda wildcards: config["models"][wildcards.model]["funcname"],
        model_params = lambda wildcards: config["models"][wildcards.model]["model_params"],
        upper_bound = lambda wildcards: config["models"][wildcards.model]["upper_bound"],
        lower_bound = lambda wildcards: config["models"][wildcards.model]["lower_bound"],
    threads: 1
    run:
        func = getattr(dadi_models, params.funcname) # extract function defining model
        #func_ex = dadi.Numerics.make_extrap_log_func(func) # wrapper to do grid extrapolation for dadi approximation of expected sfs given model
        func_ex = dadi.Numerics.make_extrap_func(func)
        
        # load data
        fs = dadi.Spectrum.from_file(input.dadisfs)
        fs = fs.fold()
        ns = fs.sample_sizes
        
        # pts is list of 3 grid sizes that are extrapolated when approximating expected sfs in the population frquency grid space (phi). dadi manual recommends setting minimum grid size to "slightly" larger than largest population sample size. But does not specify if means haploid or diploid sample size and I haven't so far tried to figure out. Seems ok to set it larger than largest haploid sample size, which for pigs is 30.
        pts = [30,40,50]

        
        pars = params.model_params
        npars = len(pars)
        lower = np.array(params.lower_bound, dtype=np.float64)
        upper = np.array(params.upper_bound, dtype=np.float64)

        p0 = np.ones(npars)
        np.random.seed(int(params.seed))
        p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper, lower_bound=lower)
        
        # this does the inference
        ml_pars = dadi.Inference.optimize_log(p0, fs, func_ex, pts,
                                   lower_bound=lower,
                                   upper_bound=upper,
                                   verbose=0, maxiter=30)

        # save maximum likelihood parameter estiamtes
        np.savetxt(output.maxlpars, ml_pars)

        # evaluate model likelihood at estimated ml params
        model = func_ex(ml_pars, ns, pts)
        like = np.array([dadi.Inference.ll_multinom(model, fs)])

        # save likelihood
        np.savetxt(output.like, like)





rule collect_fit_runs:
    input:
        runs = expand(os.path.join(OUTMAIN, "fit", "{{pair}}", "{{model}}", "run_{i}.{f}"), i=range(1, int(config["nruns"])+1), f=["maxlpars", "likehood"])
    output:
        f = os.path.join(OUTMAIN, "fit", "{pair}", "{model}", "collected_first_round_pars.txt")
    params:
        inprefix = os.path.join(OUTMAIN, "fit", "{pair}", "{model}", "run_"),
        model_params = lambda wildcards: config["models"][wildcards.model]["model_params"],
        nruns = config["nruns"]
    run:
        parnames = params.model_params
        outfile = open(output.f,"w+")
        outfile.write("seed\t" + "\t".join(parnames) + "\tlikelihood\n")

        for i in range(1, params.nruns + 1):
            infileml = params.inprefix + str(i) + ".maxlpars"
            infilelike = params.inprefix + str(i) + ".likehood"
            ml = np.genfromtxt(infileml)
            lk = np.genfromtxt(infilelike)

            outfile.write(str(i) + "\t" + "\t".join(map(str,ml)) + "\t" + str(lk) + "\n")
            
        outfile.close()




rule converge_getmaxlike:
    input:
        f = os.path.join(OUTMAIN, "fit", "{pair}", "{model}", "collected_first_round_pars.txt"),
        dadisfs = lambda wildcards: os.path.join(OUTMAIN, "2dsfs_dadi", "{p[0]}_{p[1]}_unfolded_dadi.sfs".format(p=config["pairs"][wildcards.pair])),
    output:
        f = os.path.join(OUTMAIN, "fit",  "{pair}", "{model}", "maxlike_params.txt"),
        f2 = os.path.join(OUTMAIN, "fit", "{pair}", "{model}", "collected_pars_allruns.txt"),
    params:
        maxn = config["maxn"],
        funcname = lambda wildcards: config["models"][wildcards.model]["funcname"],
        model_params = lambda wildcards: config["models"][wildcards.model]["model_params"],
        upper_bound = lambda wildcards: config["models"][wildcards.model]["upper_bound"],
        lower_bound = lambda wildcards: config["models"][wildcards.model]["lower_bound"],
    run:
        df = pd.read_table(input.f)
        # check for convergence
        converged = np.sum(abs(df["likelihood"] - np.max(df["likelihood"]))<5) >= 3
        # if not converged, load required variables to rerun. see rule fit_model for details
        if not converged:
            func = getattr(dadi_models, params.funcname) # extract function defining model
            func_ex = dadi.Numerics.make_extrap_func(func)
            fs = dadi.Spectrum.from_file(input.dadisfs)
            fs = fs.fold()
            ns = fs.sample_sizes
            pts = [30,40,50]
            pars = params.model_params
            npars = len(pars)
            lower = np.array(params.lower_bound, dtype=np.float64)
            upper = np.array(params.upper_bound, dtype=np.float64)
            i = df.shape[0] + 1

        while not converged:
            if(i > params.maxn):
                with open(output.f, "w+") as fh:
                    fh.write("did not converge after {0} runs\n see in {1} all parameters and likelihoods in each run\n".format(params.maxn, output.f2))
                break

            p0 = np.ones(npars)
            np.random.seed(i)
            p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper, lower_bound=lower)

            ml_pars = dadi.Inference.optimize_log(p0, fs, func_ex, pts,
                        lower_bound=lower,
                        upper_bound=upper,
                        verbose=0, maxiter=30)

            model = func_ex(ml_pars, ns, pts)
            like = dadi.Inference.ll_multinom(model, fs)

            # append new optimization
            df.loc[len(df)] = [x for y in [[i], ml_pars, [like]] for x in y]
            converged = np.sum(abs(df["likelihood"] - np.max(df["likelihood"]))<5) >= 3
            i += 1
            

        if converged:
            maxl = df.loc[df["likelihood"].argmax()]
            with open(output.f, "w+") as fh:
                fh.write("\t".join(maxl.index) + "\n")
                fh.write("\t".join(map(str,maxl)) + "\n")

        df.to_csv(output.f2, index=False, header=True, sep="\t")
                                 



rule rescale_plot_mlestimates:
    input:
        dadisfs = lambda wildcards: os.path.join(OUTMAIN, "2dsfs_dadi", "{p[0]}_{p[1]}_unfolded_dadi.sfs".format(p=config["pairs"][wildcards.pair])),
        mlpars = os.path.join(OUTMAIN, "fit",  "{pair}", "{model}", "maxlike_params.txt")
    output:
        fit_plot = os.path.join(OUTMAIN, "fit", "{pair}", "{model}", "{pair}_{model}_dadifit_plot.png"),
        mlpars = os.path.join(OUTMAIN, "fit", "{pair}", "{model}", "rescaled_maxlike_params.txt")
    params:
        funcname = lambda wildcards: config["models"][wildcards.model]["funcname"],
        model_params = lambda wildcards: config["models"][wildcards.model]["model_params"],
        p1 = lambda wildcards: config["pairs"][wildcards.pair][0],
        p2 = lambda wildcards: config["pairs"][wildcards.pair][1],
        mu = config["mu"],
        g = config["g"]
    run:
        with open(input.mlpars, "r") as fh:
            converged = fh.read().split()[0] == "seed"
            
        if converged:
            func = getattr(dadi_models, params.funcname) # extract function defining model
            func_ex = dadi.Numerics.make_extrap_func(func)
            fs = dadi.Spectrum.from_file(input.dadisfs)
            fs = fs.fold()
            fs.pop_ids = [params.p1, params.p2]
            ns = fs.sample_sizes
            pts = [30,40,50]
            pars = params.model_params
            npars = len(pars)

            mlres = pd.read_table(input.mlpars)
            
            ml_pars =  list(mlres.iloc[0][1:-1])
            model = func_ex(ml_pars, ns, pts)

            dadi.Plotting.plot_2d_comp_multinom(model, fs, show=False)
            plt.savefig(output.fit_plot)
            plt.clf()

            # quantities we need for the rescaling
            mu = params.mu
            g = params.g
            theta = dadi.Inference.optimal_sfs_scaling(model, fs)
            L = fs.data.sum()

            rescaled_pars = ["Na"]
            [rescaled_pars.append(x) for x in pars]

            Na = theta / (4 * mu * L)
            rescaled_ml_pars = [Na]
            for i in range(npars):
                par = pars[i]
                ml_par = ml_pars[i]
                if par[0] == "N":
                    p = ml_par * Na
                elif par[0] == "m":
                    p = ml_par / (2*Na)
                elif par[0] == "T":
                    p = ml_par * Na * 2 * g
                else:
                    p = None
                rescaled_ml_pars.append(p)
            with open(output.mlpars, "w+") as fh:
                fh.write(" ".join(rescaled_pars) + "\n")
                fh.write(" ".join(map(str, rescaled_ml_pars)) + "\n")
            

        else:
            not_converged_string = "optimization did not converge for this model.\ngo to {} for details\n".format(input.mlpars)
            plt.text(x=0.05, y=0.5, s=not_converged_string)
            plt.savefig(output.fit_plot)
            plt.clf()
            
            with open(output.mlpars, "w+") as fh:
                fh.write(not_converged_string)





rule model_selection_aic:
    input:
        mlpars = expand(os.path.join(OUTMAIN, "fit", "{{pair}}", "{model}", "maxlike_params.txt"), model= config["models"].keys())
    output:
        aic_table = os.path.join(OUTMAIN, "fit", "{pair}", "model_selection_aic.tsv")
    params:
        pair = "{pair}",
        models = config["models"].keys(),
    run:
        res = []
        for m in params.models:
            npars = len(config["models"][m]["model_params"])
            with open(os.path.join(OUTMAIN, "fit", params.pair, m, "maxlike_params.txt"), "r+") as fh:
                converged = fh.read().split()[0] == "seed"
            if not converged:
                continue

            with open(os.path.join(OUTMAIN, "fit", params.pair, m, "maxlike_params.txt"), "r+") as fh:
                mlike = float(fh.readlines()[1].split()[-1].rstrip())

            aic = -2 * mlike + 2 * npars
            res.append([m, mlike, npars, aic])

        df = pd.DataFrame(res, columns=["model", "likelihood", "npars", "AIC"])
        df.to_csv(output.aic_table, index=False, header=True, sep="\t")
