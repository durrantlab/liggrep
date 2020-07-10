LigGrep 1.0.0
=============

usage: liggrep.py [-h] [-m MODE] [-o BABEL_EXEC] [-f FILE] [-v]
                  [--num_processors N]
                  [--job_manager {serial,multiprocessing,mpi}] [-t] [-i]
                  receptor ligands filters

LigGrep 1.0.0, a free, open-source program for identifying docked poses that
participate in user-specified receptor/ligand interactions. As input, LigGrep
accepts a protein receptor file (PDB, PDBQT), a directory containing many
docked-compound files (PDB, PDBQT, SDF), and a list of user-specified filters
(JSON). It evaluates each docked pose and outputs the names of the compounds
with poses that pass all filters.

positional arguments:
  receptor              PDBQT file containing receptor information.
  ligands               Directory of PDBQT files containing the docked-ligand
                        files.
  filters               JSON file containing filters, which are formatted as a
                        list of dictionaries.

optional arguments:
  -h, --help            show this help message and exit
  -m MODE, --mode MODE  Optional user-specified bond-order mode. If OPENBABEL,
                        LigGrep will assign bond orders using the open-babel
                        executable specified via the --babel_exec parameter.
                        If SMILES, LigGrep will use SMILES files to assign
                        bond orders. These files must be in the same directory
                        as the docked-ligand files, and they must be similarly
                        named (except with the .smi extension). If NONE,
                        LigGrep will not attempt to assign bond orders. The
                        docked ligands must either be in a format that itself
                        specifies bond orders (e.g., SDF), or the user-
                        specified filters must not depend on bond orders.
                        Default: NONE
  -o BABEL_EXEC, --babel_exec BABEL_EXEC
                        Optional path to the OpenBabel executable.
  -f FILE, --file FILE  The name of the file were LigGrep analysis should be
                        saved. Defaults to "output.txt".
  -v, --verbose         Indicate why molecules and poses are rejected
                        (standard output).
  --num_processors N, -p N
                        Number of processors to use for parallel calculations.
                        Defaults to 1.
  --job_manager {serial,multiprocessing,mpi}
                        Determine what style of multiprocessing to use:
                        serial, mpi, or multiprocessing. Serial will override
                        the num_processors flag, forcing it to be one. MPI
                        mode requires mpi4py 2.1.0 or higher. Defaults to
                        serial.
  -t, --test            Run optional tests to verify that code updates don't
                        break functionality.
  -i, --internal_test   (Optional parameter used for internal testing.)

EXAMPLES OF USE:

1. Prepare a virtual library and save all 3D models to a single SDF file in
   the present directory:

python liggrep.py ./liggrep/examples/receptors/receptor.pdb \
    "./liggrep/examples/ligands/sdf/*.sdf" ./liggrep/examples/filters.json

2. Assign bond orders to PDBQT-formatted ligands using SMI files in the same
   directory:

python liggrep.py ./liggrep/examples/receptors/receptor.pdbqt \
    "./liggrep/examples/ligands/pdb/*.pdb" ./liggrep/examples/filters.json \
    --mode SMILES

3. Or use Open Babel to assign bond orders:

python liggrep.py ./liggrep/examples/receptors/receptor.pdbqt \
    "./liggrep/examples/ligands/pdb/*.pdb" ./liggrep/examples/filters.json \
    --mode OPENBABEL --babel_exec /usr/local/bin/obabel

4. By default, LigGrep saves the names of the poses that pass all filters to
   "output.txt". You can optionally specify a different output file:

python liggrep.py ./liggrep/examples/receptors/receptor.pdb \
    "./liggrep/examples/ligands/sdf/*.sdf" ./liggrep/examples/filters.json \
    --file other_output.txt

5. You can tell LigGrep to output why each ligand/pose is accepted or
   rejected:

python liggrep.py ./liggrep/examples/receptors/receptor.pdb \
    "./liggrep/examples/ligands/sdf/*.sdf" ./liggrep/examples/filters.json \
    --verbose

6. By default, LigGrep runs in serial mode, meaning it runs on only one
   processor. You can also use multiple processors to speed the process.
   Requesting -1 processors means all processors will be used.

python liggrep.py ./liggrep/examples/receptors/receptor.pdb \
    "./liggrep/examples/ligands/sdf/*.sdf" ./liggrep/examples/filters.json \
    --job_manager multiprocessing --num_processors -1
