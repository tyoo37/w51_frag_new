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
    """
    print("Running initial dendrogram script...")
    print("Runing run_initial_dendro.sh")
    result = subprocess.run('bash run_initial_dendro.sh', shell=True, check=True)
    print("STDOUT:", result.stdout)
    print("STDERR:", result.stderr)

    for region in ['W51-E', 'W51-IRS2']:
        for band in ['B3', 'B6']:
            file_path = f'tables/{region}_{band}_initial_dendro.fits'
            assert os.path.exists(file_path), f"Error: File '{file_path}' does not exist."
    

    #compare initial dendrogram with visually selected sources
    print("Running notebooks to match initial dendrogram and visually selected sources...")
    print("Running check_adam_region.ipynb")
    run_notebook('check_adam_region.ipynb', kernel_name='base')

    for region in ['W51-E', 'W51-IRS2']:
        for band in ['B3', 'B6']:
            file_path = f'pngs/{region}_{band}_adam_selected_regions.png'
            assert os.path.exists(file_path), f"Error: File '{file_path}' does not exist."

    print("Running check_nazar_region.ipynb")
    run_notebook('check_nazar_region.ipynb', kernel_name='base')

    for region in ['W51-E', 'W51-IRS2']:
        for band in ['B3', 'B6']:
            file_path = f'pngs/{region}_{band}_nazar_selected_regions.png'
            assert os.path.exists(file_path), f"Error: File '{file_path}' does not exist."
    print("Running check_taehwa_region.ipynb")
    run_notebook('check_taehwa_region.ipynb', kernel_name='base')

    for region in ['W51-E', 'W51-IRS2']:
        for band in ['B3', 'B6']:
            file_path = f'pngs/{region}_{band}_taehwa_selected_regions.png'
            assert os.path.exists(file_path), f"Error: File '{file_path}' does not exist."

    #merge the matching sources from the three coauthors into one catalog
    print("Running merge_catalog.py to merge catalogs...")
    exec(open("merge_catalog.py").read())
  
    #check whether the final dendrogram is correct
    print("Running check_final_catalog.py to verify final dendrogram...")
    exec(open('check_final_catalog.py').read())
    """
    #match the dendrogram with the band 3 and band 6 catalogs
    print("Running matching_band3band6.ipynb to match dendrogram with catalogs...")
    run_notebook('matching_band3band6.ipynb', kernel_name='base')

    #create reg files from the catalogs to check catalogs in carta
    print("Running cat_to_reg.py to create reg files from catalogs...")
    exec(open("cat_to_reg.py").read())
  
    print("Running remove_duplicates.ipynb to remove duplicate sources from the matched catalog...")
    run_notebook('remove_duplicates.ipynb', kernel_name='base')

    #manually adjust the matched catalog
    print("Running manual_adjustment.ipynb to manually adjust the matched catalog...")
    run_notebook('manual_adjustment.ipynb', kernel_name='base')
    
    print("Running make_insignificant_catalog.py to create insignificant catalog...")
    exec(open("make_insignificant_catalog.py").read())

    print("Running merge_ambiguous.ipynb to merge ambiguous sources...")
    run_notebook('merge_ambiguous.ipynb', kernel_name='base')
    """
    run_notebook('make_cutout_low_sn_and_ambiguous.ipynb', kernel_name='base')
    run_notebook('master_catalog.ipynb', kernel_name='base')
    run_notebook('making_table.ipynb', kernel_name='base')
    """