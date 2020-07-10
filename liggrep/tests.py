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

import os
import pathlib
from .__init__ import main
import sys

# Note that ligand 11 matches json file. Be sure especially to leave
# decoy.vina_out.34.pdb in there. It matches the filters if exclude: false. So
# it's a good way to test the exclusion system.


def start():
    """Start the tests."""

    example_path = (
        str(pathlib.Path(__file__).parent.absolute()) + os.sep + "examples" + os.sep
    )
    receptor_path = example_path + "receptors" + os.sep
    ligand_path = example_path + "ligands" + os.sep
    receptor_pdb = receptor_path + "receptor.pdb"
    receptor_pdbqt = receptor_path + "receptor.pdbqt"
    ligand_pdbs = ligand_path + "pdb" + os.sep + "*.pdb"
    ligand_pdbs_no_smis = ligand_path + "pdb_no_smi" + os.sep + "*.pdb"
    ligand_pdbqts = ligand_path + "pdbqt" + os.sep + "*.pdbqt"
    ligand_sdfs = ligand_path + "sdf" + os.sep + "*.sdf"
    filters = example_path + "filters.json"
    filters_pdb_ligs = example_path + "filters_pdb_ligs.json"

    obabel_paths = [f for f in sys.argv if "obabel" in f]
    if len(obabel_paths) != 1:
        # /usr/local/bin/obabel
        print(
            "FAILED: Cannot test OPENBABEL mode because no obabel executable specified at the command line "
            + "(e.g., python liggrep.py --test /usr/local/bin/obabel)"
        )
        sys.exit(0)
    obabel_path = obabel_paths[0]

    test_error_catching(
        ligand_pdbs,
        ligand_pdbqts,
        ligand_sdfs,
        ligand_pdbs_no_smis,
        filters,
        receptor_pdb,
        example_path,
        obabel_path,
    )

    test_correct_use(
        ligand_pdbs,
        ligand_pdbqts,
        ligand_sdfs,
        receptor_pdb,
        receptor_pdbqt,
        filters,
        filters_pdb_ligs,
        example_path,
        obabel_path,
    )

    print("PASSED all tests")


def detect_right_output(params, txt):
    """Detect whether a liggrep run has produced expectecd output.

    :param params: A dictionary of the user parameters and filters.
    :type params: dict
    :param txt: The text that the output must contain.
    :type txt: str
    """

    if len(params["debug_msg"].strip().split("\n")) != 1:
        print("FAILED: More than one ligand matched (but should be only one)")
        sys.exit(0)
    detect_output(params, txt, "debug_msg")


def detect_error_output(params, txt):
    """Detect whether a liggrep run has produced an expected error.

    :param params: A dictionary of the user parameters and filters.
    :type params: dict
    :param txt: The text that the error must contain.
    :type txt: str
    """

    detect_output(params, txt, "debug_error")


def detect_output(params, txt, params_key):
    """Detect whether a liggrep run has produced expected output (msg, error,
    warn, etc.).

    :param params: A dictionary of the user parameters and filters.
    :type params: dict
    :param txt: The text that the error must contain.
    :type txt: str
    :param params_key: The name of the item to examine. Like "debug_msg".
    :type params_key: str
    """

    if txt in params[params_key].replace("\n", " "):
        print('PASSED: Detected "' + txt.strip() + '"')
    else:
        print('FAILED: Did not detect "' + txt.strip() + '"')
        sys.exit(0)


def detect_missing_file(params, filename):
    """Detect whether liggrep has identified a file as being missing.

    :param params: A dictionary of the user parameters and filters.
    :type params: dict
    :param filename: The filename.
    :type filename: str
    """

    if "ERROR" in params["debug_error"] and filename in params["debug_error"]:
        print("PASSED: Detected absent file/mask: " + filename)
    else:
        print("FAILED: Did not detect absent file/mask: " + filename)
        sys.exit(0)


def test_error_catching(
    ligand_pdbs,
    ligand_pdbqts,
    ligand_sdfs,
    ligand_pdbs_no_smis,
    filters,
    receptor_pdb,
    example_path,
    obabel_path,
):
    """Test error catching.

    :param ligand_pdbs: A mask describing the pdb ligand files.
    :type ligand_pdbs: str
    :param ligand_pdbqts: A mask describing the pdbqt ligand files.
    :type ligand_pdbqts: str
    :param ligand_sdfs: A mask describing the sdf ligand files.
    :type ligand_sdfs: str
    :param ligand_pdbs_no_smis: A mask describing the pdb ligand files without
        associated smi files.
    :type ligand_pdbs_no_smis: str
    :param filters: The path to the filters JSON file.
    :type filters: str
    :param receptor_pdb: The path to the receptor PDB file.
    :type receptor_pdb: str
    :param example_path: The path to the example files.
    :type example_path: str
    :param obabel_path: The path to the obabel executable.
    :type obabel_path: str
    """

    # Test the various file-exists errors.
    print("TEST: Missing receptor")
    p = main(["__no_exist.receptor.pdb", ligand_pdbs, filters, "-i"])
    detect_missing_file(p, "__no_exist.receptor.pdb")

    print("TEST: Missing ligand PDB files")
    p = main([receptor_pdb, "./__no_files/*.pdb", filters, "-i"])
    detect_missing_file(p, "./__no_files/*.pdb")

    print("TEST: Missing filters file")
    p = main([receptor_pdb, ligand_pdbs, "__no_exist.filters.json", "-i"])
    detect_missing_file(p, "__no_exist.filters.json")

    print('TEST: Both "receptorAtom" and "coordinate" in filter')
    p = main([receptor_pdb, ligand_pdbs, example_path + "filters_bad.json", "-i"])
    detect_error_output(p, 'One of your filters contains both a "receptorAtom"')

    print("TEST: All filter keys are valid")
    p = main([receptor_pdb, ligand_pdbs, example_path + "filters_bad3.json", "-i"])
    detect_error_output(p, "contains an invalid key")

    print("TEST: All filters contain required keys")
    p = main([receptor_pdb, ligand_pdbs, example_path + "filters_bad4.json", "-i"])
    detect_error_output(p, "does not contain the required key")

    print("TEST: All filters have query-point specification")
    p = main([receptor_pdb, ligand_pdbs, example_path + "filters_bad5.json", "-i"])
    detect_error_output(p, "but not both")

    print("TEST: Receptor atom does not exist")
    p = main([receptor_pdb, ligand_pdbs, example_path + "filters_bad2.json", "-i"])
    detect_error_output(p, "No receptor atom matches filter")

    print("TEST: Bad user-specified SMARTS string")
    p = main([receptor_pdb, ligand_sdfs, example_path + "filters_bad6.json", "-i"])
    detect_error_output(p, "SMARTS string is not valid: [#8](sd323[")

    for mode in ["OPENBABEL", "SMILES"]:
        print("TEST: Test error if SDF file with " + mode + " mode.")
        p = main(
            [
                receptor_pdb,
                ligand_sdfs,
                filters,
                "-m",
                "SMILES",
                "--babel_exec",
                obabel_path,
                "-i",
            ]
        )
        detect_error_output(p, "OPENBABEL/SMILES mode is not appropriate for SDF files")

    for ligs, ligs_desc in [
        (ligand_pdbs, "PDB"),
        (ligand_pdbqts, "PDBQT"),
    ]:
        print(
            (
                "TEST: Trying to use {ligs_desc} ligands in NONE mode with "
                + "double-bonds in filters."
            ).format(ligs_desc=ligs_desc)
        )
        p = main([receptor_pdb, ligs, filters, "-m", "NONE", "-i",])
        detect_error_output(p, "PDB- and PDBQT-formatted ligands in NONE")

    print("TEST: OPENBABEL specified, but no obabel path given")
    p = main([receptor_pdb, ligand_pdbs, filters, "-m", "OPENBABEL", "-i",])
    detect_error_output(p, "via the --babel_exec parameter")

    print("TEST: Open babel path does not exist")
    p = main(
        [
            receptor_pdb,
            ligand_pdbs,
            filters,
            "-m",
            "OPENBABEL",
            "-i",
            "--babel_exec",
            "/fake/path/obabel",
        ]
    )
    detect_error_output(p, '--babel_exec parameter does not exist: "/fake/path/obabel"')

    print("TEST: Open babel path is not obabel (but perhaps babel)")
    p = main(
        [
            receptor_pdb,
            ligand_pdbs,
            filters,
            "-m",
            "OPENBABEL",
            "-i",
            "--babel_exec",
            "/fake/path/babel",
        ]
    )
    detect_error_output(p, 'paramter does not start with "obabel"')

    print("TEST: No SMI files present when running in SMILES mode.")
    p = main([receptor_pdb, ligand_pdbs_no_smis, filters, "-m", "SMILES", "-i"])
    detect_error_output(p, "not have an associated SMILES")


def test_correct_use(
    ligand_pdbs,
    ligand_pdbqts,
    ligand_sdfs,
    receptor_pdb,
    receptor_pdbqt,
    filters,
    filters_pdb_ligs,
    example_path,
    obabel_path,
):
    """Additional tests.

    :param ligand_pdbs: A mask describing the pdb ligand files.
    :type ligand_pdbs: str
    :param ligand_pdbqts: A mask describing the pdbqt ligand files.
    :type ligand_pdbqts: str
    :param ligand_sdfs: A mask describing the sdf ligand files.
    :type ligand_sdfs: str
    :param receptor_pdb: The path to the receptor PDB file.
    :type receptor_pdb: str
    :param receptor_pdbqt: The path to the receptor PDBQT file.
    :type receptor_pdbqt: str
    :param filters: The path to the filters JSON file.
    :type filters: str
    :param filters_pdb_ligs: The path to the filters JSON file appropriate for
        PDB ligand files.
    :type filters_pdb_ligs: str
    :param example_path: The path to the example files.
    :type example_path: str
    :param obabel_path: The path to the obabel executable.
    :type obabel_path: str
    """

    # Test NONE mode
    for recep, recep_desc in [(receptor_pdb, "PDB"), (receptor_pdbqt, "PDBQT")]:
        for ligs, constrnt_file, ligs_desc in [
            (ligand_pdbs, filters_pdb_ligs, "PDB"),
            (ligand_sdfs, filters, "SDF"),
            (ligand_pdbqts, filters_pdb_ligs, "PDBQT"),
        ]:
            print(
                "TEST: {recep_desc} receptor + {ligs_desc} ligand (NONE mode)".format(
                    recep_desc=recep_desc, ligs_desc=ligs_desc
                )
            )
            # print(" ".join([recep, ligs, constrnt_file, "-m", "NONE", "-i"]))
            p = main([recep, ligs, constrnt_file, "-m", "NONE", "-i"])
            detect_right_output(
                p, "decoy.vina_out.11." + ligs_desc.lower() + " (pose 3)",
            )

    # Test OPENBABEL mode (pdb and pdbqt)
    for recep, recep_desc in [(receptor_pdb, "PDB"), (receptor_pdbqt, "PDBQT")]:
        for ligs, constrnt_file, ligs_desc in [
            (ligand_pdbs, filters, "PDB"),
            (ligand_pdbqts, filters, "PDBQT"),
        ]:
            print(
                (
                    "TEST: {recep_desc} receptor + {ligs_desc} ligand "
                    + "(OPENBABEL mode)"
                ).format(recep_desc=recep_desc, ligs_desc=ligs_desc)
            )
            p = main(
                [
                    recep,
                    ligs,
                    constrnt_file,
                    "-m",
                    "OPENBABEL",
                    "--babel_exec",
                    obabel_path,
                    "-i",
                ]
            )
            detect_right_output(
                p, "decoy.vina_out.11." + ligs_desc.lower() + " (pose 3)",
            )

    # Test SMILES mode
    for recep, recep_desc in [(receptor_pdb, "PDB"), (receptor_pdbqt, "PDBQT")]:
        for ligs, constrnt_file, ligs_desc in [
            (ligand_pdbs, filters, "PDB"),
            (ligand_pdbqts, filters, "PDBQT"),
        ]:
            print(
                (
                    "TEST: {recep_desc} receptor + {ligs_desc} ligand "
                    + "(SMILES mode)"
                ).format(recep_desc=recep_desc, ligs_desc=ligs_desc)
            )
            # print(" ".join([recep, ligs, constrnt_file, "-m", "SMILES", "-i"]))
            p = main([recep, ligs, constrnt_file, "-m", "SMILES", "-i"])
            detect_right_output(
                p, "decoy.vina_out.11." + ligs_desc.lower() + " (pose 3)",
            )

    # Make sure it still works when you only specify partial receptor ligand
    # atoms.
    print(
        "TEST: PDB receptor + PDB ligand "
        + "(OPENBABEL mode) with partial receptor atoms specified"
    )
    p = main(
        [
            receptor_pdb,
            ligand_pdbs,
            example_path + "filters_limit_selection.json",
            "-m",
            "OPENBABEL",
            "--babel_exec",
            obabel_path,
            "-i",
        ]
    )
    detect_right_output(
        p, "decoy.vina_out.11.pdb (pose 3)",
    )

    # Test OPENBABEL mode on pdb ligands using all processors (to test
    # multiprocessing).
    print(
        "TEST: PDB receptor + PDB ligand "
        + "(OPENBABEL mode) using multiple (all) processors "
        + "(final test, will fail on Windows)"
    )
    p = main(
        [
            receptor_pdb,
            ligand_pdbs,
            filters,
            "-m",
            "OPENBABEL",
            "--babel_exec",
            obabel_path,
            "-i",
            "--num_processors",
            "-1",
            "--job_manager",
            "multiprocessing",
        ]
    )
    detect_right_output(
        p, "decoy.vina_out.11.pdb (pose 3)",
    )
