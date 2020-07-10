# Copyright 2020 Jacob D. Durrant
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not
# use this file except in compliance with the License. You may obtain a copy
# of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
# WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
# License for the specific language governing permissions and limitations
# under the License.

from subprocess import Popen, PIPE
import os
from .. import output
from . import load


def get_ligand_mols(ligand_filename, babel_exec, params):
    """
       Create a molecule of the ligand, with the appropriate hybridization
       states assigned, using Open Babel.

       :param ligand_filename: File name of the ligand pdb or pdbqt file.
       :type ligand_filename: str
       :param babel_exec: Path to Open Babel executable.
       :type babel_exec: str
       :param params: A dictionary of the user parameters and filters.
       :type params: dict
       :returns: A list of RDKit molecule objects, containing the
           appropriate atom information and hybridization states of the ligand
           conformations.
       :rtype: list
    """

    filename_root, filename_ext = os.path.splitext(ligand_filename)

    ligand_mols = None

    # First, try to convert using open babel and the -h flag, then -p if that
    # fails.
    for flag in ["-h", "-p"]:
        obabel_convert(filename_root, filename_ext, babel_exec, flag)
        ligand_mols = load.load_ligand_mols(filename_root + ".tmp.sdf", params)

        # Delete SDF.
        os.unlink(filename_root + ".tmp.sdf")

        if ligand_mols is not None:
            break

    if ligand_mols is None:
        output.warn("Could not convert " + ligand_filename + " to SDF.", params)

    return ligand_mols


def obabel_convert(filename_root, filename_ext, babel_exec, flag):
    """
       Convert the ligand pdb or pdbqt file to a sdf file using Open Babel.
       This captures the appropriate hybridization states of the ligand molecule.

       :param filename_root: Root name of the ligand pdb or pdbqt file.
       :type filename_root: str
       :param filename_ext: Extension of the ligand file, either pdb or pdbqt.
       :type filename_ext: str
       :param babel_exec: Path to Open Babel executable.
       :type babel_exec: str
       :param flag: The obabel -h parameter adds hydrogens, while the obabel -p
           (pH) parameter ionizes molecules as appropriate for a default pH of 7.4.
       :type flag: str
    """

    output_flnm = filename_root + ".tmp.sdf"

    if os.path.exists(output_flnm):
        os.unlink(output_flnm)

    if filename_ext == ".pdbqt":
        process = Popen(
            [
                babel_exec,
                filename_root + ".pdbqt",
                "-O",
                output_flnm,
                flag,
                " ---errorlevel 1",  # Note that default is 2
            ],
            stdout=PIPE,
            stderr=PIPE,
        )
        stdout, stderr = process.communicate()
        # if (stderr.strip() != ""):
        #     print(stderr)
        process.wait()
    elif filename_ext == ".pdb":
        process = Popen(
            [
                babel_exec,
                filename_root + ".pdb",
                "-O",
                output_flnm,
                flag,
                " ---errorlevel 1",  # Note that default is 2
            ],
            stdout=PIPE,
            stderr=PIPE,
        )
        stdout, stderr = process.communicate()
        # if (stderr.strip() != ""):
        #     print(stderr)
        process.wait()
