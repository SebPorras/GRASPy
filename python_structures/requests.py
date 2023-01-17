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
from file_formatter import *


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

    params["Tree"] = nwkToJSON(tree)

    params["Alignment"] = readAln(aln, "Protein")

    request["Params"] = params

    j_request = json.dumps(request) + '\n'

    return client.sendRequest(j_request)


def requestJointReconstruction(aln: str, nwk: str,
                               auth: str = "Guest",
                               indels: str = "BEP",
                               model: str = "JTT"):
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

    params["Tree"] = nwkToJSON(tree)
    params["Alignment"] = readAln(aln, "Protein")

    params["Inference"] = "Joint"
    params["Indels"] = indels
    params["Model"] = model

    request["Params"] = params

    j_request = json.dumps(request) + '\n'

    return client.sendRequest(j_request)


tree = requestJointReconstruction(aln="./test_data/big_test_data/GRASPTutorial_Final.aln",
                                  nwk="./test_data/big_test_data/GRASPTutorial_Final.nwk")


# tree = requestJointReconstruction(aln="./test_data/big_test_data/GRASPTutorial_Final.aln",
#                                   nwk="./test_data/big_test_data/_ancestors.nwk")

# print(tree)
# requestPlaceInQueue(29)

# output = requestJobOutput(32)
