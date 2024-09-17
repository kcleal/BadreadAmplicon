"""
Simulate fusion-seq reads with varying amplification (Nanopore long-read sequencing)

Inputs:
-------
folder of fasta files created by generate_fusions2

Outputs:
--------
fastg files of sequenced telomere fusions


"""

import subprocess
import click
from multiprocessing import Pool
import os

@click.command()
@click.option("--reference", help="path to the folder containing generated fusions", required=True)
@click.option("--out", help="output folder", required=True)
@click.option("--threads", default=1, help="number of cores", type=int)
def simulate_clusters(**args):
    ref_dir = args['reference']
    out_dir = args['out']

    files_to_process = [f for f in os.listdir(ref_dir) if 'fusions' in f]

    arguments = [(os.path.join(ref_dir, file), out_dir) for file in files_to_process]  # Full path for files
    threads = args['threads']
    with Pool(threads) as pool:
        pool.starmap(run_simulation, arguments)


def run_simulation(file_path, out_dir):
    filename = os.path.basename(file_path)
    quantity = filename.split('_')[1]
    name = filename.split('.')[0]
    output_file = os.path.join(out_dir, f'{name}.fastq')
    command = ['badread', 'simulate', '--reference', file_path, '--quantity', f'{quantity}x']

    try:
        with open(output_file, 'wb') as out_f:
            result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
            out_f.write(result.stdout)
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {' '.join(command)}")
        if e.stderr:
            print(f"Error message: {e.stderr.decode()}")
    except Exception as e:
        print(f"Unexpected error: {str(e)}")


if __name__ == "__main__":
    simulate_clusters()

