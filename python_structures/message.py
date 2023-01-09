############################################
# Date: 9/1/23
# Author: Sebastian Porras
# Aims: Message class to convert aln and nwk to JSON for socket
#####################################################################

import json


def serialiseAln(file_name: str) -> dict[str, str]:
    """
    Creates a dictionary where seq ids are the key 
    and alignment is the value. 

    Parameters:
        file_name(str): path to aln file 

    Returns:
        dict   type="Joint"
)

    """

    alns = {}

    with open(file_name, "r") as fa:

        line = fa.readline()

        while line:

            if line[0] == ">":

                key = "_".join(line.split())

                aln = ""

                line = fa.readline()

                while line[0] != ">":

                    aln += line.strip()
                    line = fa.readline()

                    if not line:
                        break

                alns[key] = aln

            else:
                line = fa.readline()

    return alns


def constructJsonMessage(aln: str, nwk: str, type: str) -> str:

    message = dict()

    with open(nwk, 'r') as f:
        tree = ""
        for line in f:
            tree += line.strip()

    message["tree"] = tree
    message["alignments"] = serialiseAln(aln)
    message["type"] = type

    return json.dumps(message)
