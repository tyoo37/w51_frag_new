

"""
https://stackoverflow.com/questions/67811531/how-can-i-execute-a-ipynb-notebook-file-in-a-python-script
https://nbconvert.readthedocs.io/en/latest/execute_api.html
https://github.com/keflavich/brick-jwst-2221/blob/3345a8f7d2bdab9b09dac29dbbff67632dd0ec5d/brick2221/reduction/run_notebook.py
"""
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor
import subprocess
import os

def run_notebook(filename, kernel_name='python39'):
    print(f"Running notebook {filename}")
    with open(filename) as ff:
        nb_in = nbformat.read(ff, nbformat.NO_CONVERT)

    print(f"Writing backup notebook {filename}.backup")
    with open(filename + ".backup", 'w', encoding='utf-8') as fh:
        nbformat.write(nb_in, fh)

    ep = ExecutePreprocessor(timeout=1200)

    nb_out = ep.preprocess(nb_in)

    print(f"Writing notebook {filename}")
    with open(filename, 'w', encoding='utf-8') as fh:
        nbformat.write(nb_in, fh)

    return nb_in, nb_out

if __name__ == "__main__":
    #run dendrogram
    os.chdir('/home/t.yoo/w51/w51_frag_new/dendro')
    print("Running dendrogram all process script...")
    #exec(open("dendro/dendro_all_process_script.py").read())

    # run convolve script
    print("Running convolve script...")
    #exec(open("convolve/convolve.py").read())

    #run TGIF to create fitting results
    os.chdir('/home/t.yoo/w51/w51_frag_new/flux')

    print("Running TGIF script...")
    #exec(open("run_tgif.py").read())    
    
    os.chdir('/home/t.yoo/w51/w51_frag_new/size')
    #get size distribution 
    print("Running size distribution notebook...")
    #run_notebook('size_distribution.ipynb', kernel_name='base')

    #run spectral index notebook / will opbtain spetral index
    os.chdir('/home/t.yoo/w51/w51_frag_new/spectral_index')
    print("Running spectral index notebook...")
    #run_notebook('spectral_index.ipynb', kernel_name='base')

    # will obtain temperature
    print("Running flux flux notebook...")
    #run_notebook('flux_flux_new.ipynb', kernel_name='base')

    # will obtain flux mass histogram and use temperature and flux
    os.chdir('/home/t.yoo/w51/w51_frag_new/flux')
    print("Running flux mass histogram notebook...")
    #run_notebook('flux_mass_histogram.ipynb', kernel_name='base')
    print("Running main fragmentation notebook...")
    # this notebook needs to be refactored in the future
    os.chdir('/home/t.yoo/w51/w51_frag_new/fragment')
    run_notebook('multiplicity_newnewnewnew_fabien_manual.ipynb', kernel_name='base')
    print("Running ciPPOs vs caPPOs notebook...")
    run_notebook('isolated_vs_core_associated_copy.ipynb', kernel_name='base')


        