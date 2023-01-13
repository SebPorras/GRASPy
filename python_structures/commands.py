############################################
# Date: 9/1/23
# Author: Sebastian Porras
# Aims: params class to convert aln and nwk to JSON for socket
#####################################################################

############################################
# Date: 10/1/23
# Author: Sebastian Porras
# Aims: Starting to create formats for requests to the bnkit server
#####################################################################

import json
import client

###### FORMATTING######


def serialiseAln(file_name: str) -> dict:
    """
    Creates a dictionary where seq ids are the key
    and alignment is the value.

    Parameters:
        file_name(str): path to aln file

    Returns:
        list: alignment seqs in JSON format

    """

    seqs = []

    with open(file_name, "r") as fa:

        line = fa.readline()

        while line:

            if line[0] == ">":

                tmp = {}

                # this would save extra info in fasta file
                # key = "_".join(line.split())

                # currently, just saves seq id and removes '>'
                key = line.split()[0][1:]

                aln = ""

                line = fa.readline()

                while line[0] != ">":

                    aln += line.strip()
                    line = fa.readline()

                    if not line:
                        break

                tmp["name"] = key
                tmp["Seq"] = aln

                seqs.append(tmp)
            else:
                line = fa.readline()

    alignments = dict()
    alignments["Sequences"] = seqs
    alignments["Datatype"] = {"Predef": "Protein"}

    return alignments


###### REQUESTS######

def requestJobOutput(job_id: int) -> str:
    """Requests the output of a submitted job. Request will be
    denied if the job is not complete.

    Parameters:
        job_id(int): The ID of the job

    Returns:
        str: {"Job":<job-number>, "Result":{<result-JSON>}}
    """

    request = dict()

    request["Command"] = "Output"
    request["Job"] = job_id

    j_request = json.dumps(request) + '\n'

    return client.sendRequest(j_request)


def requestPlaceInQueue(job_id: int) -> str:
    """Requests the status of a submitted job

    Parameters:
        job_id(str): The ID of the job

    Returns:
        str: bnkit job status
    """

    request = dict()

    request["Command"] = "Place"
    request["Job"] = job_id

    j_request = json.dumps(request) + '\n'

    return client.sendRequest(j_request)


def requestCancelJob(job_id: int) -> str:
    """Requests the status of a submitted job

    Parameters:
        job_id(str): The ID of the job

    Returns:
        str: bnkit job status
    """

    request = dict()

    request["Command"] = "Retrieve"
    request["Job"] = job_id

    j_request = json.dumps(request) + '\n'

    return client.sendRequest(j_request)


###### COMMANDS######


def requestPOGTree(aln: str, nwk: str, auth: str = "Guest") -> str:
    """Queries the server to turn an alignment 
    and a nwk string of a tree into the POGTree format.

    This output JSON can be converted into a POGTree object 
    via POGTreeFromJSON()

    Parameters:
        aln(str) = file name of aln file
        nwk(str) = file name of nwk file

    Returns:
        str: JSON string of POGTree 
    """

    request = dict()

    request["Command"] = "Pogit"
    request["Auth"] = auth

    params = dict()

    with open(nwk, 'r') as f:
        tree = ""
        for line in f:
            tree += line.strip()

    params["Tree"] = tree

    params["Alignment"] = serialiseAln(aln)

    request["Params"] = params

    j_request = json.dumps(request) + '\n'

    return client.sendRequest(j_request)


def requestJointReconstruction(aln: str, nwk: str,
                               auth: str = "Guest",
                               indels: str = "BEP",
                               model: str = "JTT") -> str:
    """Queries the bnkit server for a joint reconstruction.
    Will default to standard bnkit reconstruction parameters which
    use BEP for indels and JTT for the substitution model.

    Parameters:
        aln(str) = file name of aln file
        nwk(str) = file name of nwk file
        auth(str) = Authentication token, defaults to Guest
        indels(str) = Indel mode, defaults to BEP
        model(str) = Substitution model, defaults to JTT

    Returns:
        str: {"Message":"Queued","Job":<job-number>}
    """

    request = dict()

    request["Command"] = "Recon"
    request["Auth"] = auth

    params = dict()

    with open(nwk, 'r') as f:
        tree = ""
        for line in f:
            tree += line.strip()

    params["Tree"] = tree
    params["Alignment"] = serialiseAln(aln)

    params["Inference"] = "Joint"
    params["Indels"] = indels
    params["Model"] = model

    request["Params"] = params

    j_request = json.dumps(request) + '\n'

    return client.sendRequest(j_request)


# joint = requestJointReconstruction(
#     aln="./test_data/small_test_data/test_aln.aln",
#     nwk="./test_data/small_test_data/test_nwk.nwk")

# params["Tree"] = {"Parents": [-1, 0, 1, 1, 3, 4, 4, 6, 6, 3, 9, 10, 10, 9, 0], "Labels": ["0", "1", "A", "2", "3", "B", "4", "C", "D","5", "6", "E", "F", "G", "H"], "Distances": [0, 0.5, 0.6, 3.2, 5, 3.3, 1.8, 1, 2.5, 7, 2.5, 3.9, 4.5, 0.3, 1.1], "Branchpoints": 15}


tree = requestPOGTree(aln="./test_data/small_test_data/test_aln.aln",
                      nwk="./test_data/small_test_data/test_nwk.nwk")


print(tree)
