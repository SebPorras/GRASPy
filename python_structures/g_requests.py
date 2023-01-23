
###############################################################################
# Date: 20/1/23
# Author: Sebastian Porras
# Aims: This file contains all the protocols that the client can use to
# obtain information about a job or the server. It only contains the requests
# that the user can ask of the server.
###############################################################################

import json
import client
from parsers import *


###### REQUESTS######


def requestJobOutput(job_id: int) -> dict:
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

    response = client.sendRequest(j_request)

    return json.loads(response)


def requestPlaceInQueue(job_id: int) -> str:
    """Requests the status of a submitted job

    Parameters:
        job_id(str): The ID of the job

    Returns:
        str: {"Job":<job-number>, "Place":{<place>}}
    """

    request = dict()

    request["Command"] = "Place"
    request["Job"] = job_id

    j_request = json.dumps(request) + '\n'

    response = client.sendRequest(j_request)

    return response


def requestCancelJob(job_id: int) -> str:
    """Requests the status of a submitted job

    Parameters:
        job_id(str): The ID of the job

    Returns:
        str: {"Job":<job-number>}
    """

    request = dict()

    request["Command"] = "Retrieve"
    request["Job"] = job_id

    j_request = json.dumps(request) + '\n'

    return client.sendRequest(j_request)


def requestViewQueue() -> dict:
    """Lists all the jobs currently being 
    performed by the server 


    Returns:
        str: summary of current jobs in the server 
    """

    request = dict()

    request["Command"] = "Status"

    j_request = json.dumps(request) + '\n'

    response = client.sendRequest(j_request)

    return json.loads(response)

###### COMMANDS######


def requestPOGTree(aln: str, nwk: str, auth: str = "Guest") -> dict:
    """Queries the server to turn an alignment
    and a nwk file into the POGTree format with POGraphs for extants.

    This output JSON can be converted into a POGTree object
    via POGTreeFromJSON()

    Parameters:
        aln(str) = path to file name of aln 
        nwk(str) = path to or file name of nwk 

    Returns:
        dict: Will complete the job and provide a POG graph of the
        extants and a tree or will provide the job number if queued. 
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

    params["Alignment"] = alnToJSON(aln, "Protein")

    request["Params"] = params

    j_request = json.dumps(request) + '\n'

    j_response = client.sendRequest(j_request)

    response = json.loads(j_response)

    return response


def JointReconstruction(aln: str, nwk: str,
                        auth: str = "Guest",
                        indels: str = "BEP",
                        model: str = "JTT",
                        alphabet: str = "Protein") -> str:
    """Queries the bnkit server for a joint reconstruction.
    Will default to standard bnkit reconstruction parameters which
    use BEP for indels and JTT for the substitution model.

    Currently ONLY implemented for protein alphabets. Future version
    will guess unless user specifies alphabet. 

    Parameters:
        aln(str) = path to file name of aln 
        nwk(str) = path to file name of nwk 
        auth(str) = Authentication token, defaults to Guest
        indels(str) = Indel mode, defaults to BEP
        model(str) = Substitution model, defaults to JTT
        alphabet(str) = Sequence type. e.g. DNA or Protein 

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
    params["Alignment"] = alnToJSON(aln, alphabet)

    params["Inference"] = "Joint"
    params["Indels"] = indels
    params["Model"] = model

    request["Params"] = params

    j_request = json.dumps(request) + '\n'

    response = client.sendRequest(j_request)

    return response
