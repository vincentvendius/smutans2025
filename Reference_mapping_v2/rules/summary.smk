from itertools import chain
from functools import reduce
import pathlib
import os

def is_non_zero_file(fpath):
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0


def hetsum_inputs(wildcards):
    files = expand(
        config["pdir"] + "/tables/" + config["prefix"] + "/{smp}.bestmap.het.tsv",
        smp=SAMPLENAMES,
    )
    return files



def fast_flatten(input_list):
    return list(chain.from_iterable(input_list))


def remove_suffix(input_string, suffix):
    if suffix and input_string.endswith(suffix):
        return input_string[: -len(suffix)]
    return input_string



def read_tp_stats_tables(file, suffix):
    if is_non_zero_file(file):
        #smp = remove_suffix(pathlib.Path(file).stem, suffix)
        df = pd.read_csv(file, sep="\t")
        #df.rename(columns = {smp:sample_label_dict_read[smp]}, inplace = True)
        return df


rule summary:
    input:
        hetsum_files=lambda wc: hetsum_inputs(wc),
    output:
        topmap=config["pdir"] + "/summary/" + config["prefix"] + "/topmap.tsv" 
    params:
        tablepath=config["pdir"] + "/tables/" + config["prefix"],
        fastapath=config["pdir"] + "/fasta/" + config["prefix"]
    threads: 1
    message:
        """--- Summarize taxonomic profiling."""
    run:
        df = pd.concat( #create summary stat table across all samples
            map(
                lambda file: read_tp_stats_tables(file, suffix=".bestmap.het.tsv"),
                fast_flatten(list({input.hetsum_files})),
            )
        )
        df.to_csv( #write output
            output.topmap,
            sep="\n",
            #compression="gzip",
            header=True,
            index=False,
        )
"""
        tablepath=params.tablepath
        files=os.listdir(tablepath) #get list of files in directory 
        fastapath=params.fastapath
        fastafiles=os.listdir(fastapath)
        print(files)
        for hetfile in fast_flatten(list({input.hetsum_files})): #for each sample
            print(hetfile)
            if is_non_zero_file(hetfile): #get hash number for sample 
                smp = remove_suffix(pathlib.Path(hetfile).stem,".bestmap.het")

            #Change hash name to original name in heterogenity file
            f=open(os.path.join(tablepath,smp+".bestmap.genes.tsv"),"r")
            filedata=f.read()
            f.close()
            orgname=filedata.replace(smp,sample_label_dict_read[smp])
            f=open(os.path.join(tablepath,smp+".bestmap.genes.tsv"),"w")
            f.write(orgname)
            f.close()
            #change filenames in tables dir
            for tablefile in files:
                fsuffix="."+'.'.join(tablefile.split(".")[1:]) 
                print(fsuffix)
                print(os.path.join(tablepath, tablefile))
                print(sample_label_dict_read[smp])
                os.rename(os.path.join(tablepath, tablefile), os.path.join(tablepath, ''.join([sample_label_dict_read[smp], fsuffix])))
                print("it just work")

            for fastafile in fastafiles: 
                fsuffix="."+'.'.join(fastafile.split(".")[1:]) 
                os.rename(os.path.join(fastapath, fastafile), os.path.join(fastapath, ''.join([sample_label_dict_read[smp], fsuffix])))
"""
