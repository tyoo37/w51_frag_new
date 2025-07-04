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
    #run shell script to run initial dendrogram
    print("Running initial dendrogram script...")
    print("Runing run_initial_dendro.sh")
    result = subprocess.run('run_initial_dendro.sh', shell=True, check=True)
    print("STDOUT:", result.stdout)
    print("STDERR:", result.stderr)

    assert os.path.exists(file_path), f"Error: File '{file_path}' does not exist."
    assert os.path.exists(file_path), f"Error: File '{file_path}' does not exist."
    assert os.path.exists(file_path), f"Error: File '{file_path}' does not exist."
    assert os.path.exists(file_path), f"Error: File '{file_path}' does not exist."


    #run notebook
    print("Running notebooks to match initial dendrogram and visually selected sources...")
    print("Running check_adam_regions.ipynb")
    run_notebook('check_adam_regions.ipynb', kernel_name='base')
    print("Running check_nazar_regions.ipynb")
    run_notebook('check_nazar_regions.ipynb', kernel_name='base')
    print("Running check_taehwa_regions.ipynb")
    run_notebook('check_taehwa_regions.ipynb', kernel_name='base')

    

