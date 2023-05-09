
###############################################################################
# Date: 20/1/23
# Author: Sebastian Porras
# Aims: This file contains all the protocols that the client can use to
# obtain information about a job or the server. It only contains the requests
# that the user can ask of the server.
###############################################################################

import json
from . import client
from . import parsers
from typing import Optional


###### REQUESTS######

def send_and_recieve(request: dict) -> dict:

    j_request = json.dumps(request) + '\n'

    j_response = client.sendRequest(j_request)

    response = json.loads(j_response)

    print(response)

    return response


def JobOutput(job_id: int) -> dict:
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

    # won't use function so entire output is not printed
    j_request = json.dumps(request) + '\n'

    j_response = client.sendRequest(j_request)

    response = json.loads(j_response)

    return response


def PlaceInQueue(job_id: int) -> dict[str, int]:
    """Requests the status of a submitted job

    Parameters:
        job_id(str): The ID of the job

    Returns:
        str: {"Job":<job-number>, "Place":{<place>}}
    """

    request = dict()

    request["Command"] = "Place"
    request["Job"] = job_id

    return send_and_recieve(request)


def CancelJob(job_id: int) -> dict[str, int]:
    """Requests the status of a submitted job

    Parameters:
        job_id(str): The ID of the job

    Returns:
        str: {"Job":<job-number>}
    """

    request = dict()

    request["Command"] = "Retrieve"
    request["Job"] = job_id

    return send_and_recieve(request)


def ViewQueue() -> dict:
    """Lists all the jobs currently being 
    performed by the server 


    Returns:
        str: summary of current jobs in the server 
    """

    request = dict()

    request["Command"] = "Status"

    return send_and_recieve(request)


def JobStatus(job_id: int) -> dict:
    """Retrives job status 


    Returns:
        str: status of job as completed or queued
    """

    request = dict()

    request["Command"] = "Status"
    request["Job"] = job_id

    return send_and_recieve(request)

###### COMMANDS######


def ExtantPOGTree(aln: str, nwk: str, auth: str = "Guest") -> dict:
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

    params["Tree"] = parsers.nwkToJSON(tree)

    params["Alignment"] = parsers.alnToJSON(aln, "Protein")

    request["Params"] = params

    return send_and_recieve(request)


def JointReconstruction(aln: str, nwk: str,
                        auth: str = "Guest",
                        indels: str = "BEP",
                        model: str = "JTT",
                        alphabet: Optional[str] = None) -> dict:
    """Queries the bnkit server for a joint reconstruction.
    Will default to standard bnkit reconstruction parameters which
    use BEP for indels and JTT for the substitution model.

    Current accepted alphabets: 'DNA', 'RNA', 'Protein'

    Parameters:
        aln(str) = path to file name of aln 
        nwk(str) = path to file name of nwk 
        auth(str) = Authentication token, defaults to Guest
        indels(str) = Indel mode, defaults to BEP
        model(str) = Substitution model, defaults to JTT
        alphabet(str) = Sequence type. e.g. DNA or Protein. 
                        If user does not specify, it will guess
                        based on sequence content. 

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

    params["Tree"] = parsers.nwkToJSON(tree)
    params["Alignment"] = parsers.alnToJSON(aln, alphabet)

    params["Inference"] = "Joint"
    params["Indels"] = indels
    params["Model"] = model

    request["Params"] = params

    return send_and_recieve(request)


def LearnLatentDistributions(nwk: str,
                             states: list[str],
                             csv_data: str,
                             auth: str = "Guest"
                             ) -> dict:
    """Learns the distribution of an arbitrary number of discrete 
    states. The output from the job will be a new/refined distribution. 

    ONLY implemented for continuous values not discrete multinomial 
    symbols.

    Parameters:
        nwk(str) = path to file name of nwk 
        states(list) = a list of names for each latent states
        csv_data(str) = path to csv with data 
        auth(str) = Authentication token, defaults to Guest

    Returns:
        str: {"Message":"Queued","Job":<job-number>}
    """

    request = dict()

    request["Command"] = "Train"
    request["Auth"] = auth

    params = dict()

    params["States"] = states

    # format tree
    with open(nwk, 'r') as f:
        tree = ""
        for line in f:
            tree += line.strip()

    params["Tree"] = parsers.nwkToJSON(tree)

    params["Dataset"] = parsers.csvDataToJSON(csv_data)

    # load all parameters
    request["Params"] = params

    return send_and_recieve(request)


def MarginaliseDistOnAncestor(nwk: str,
                              states: list[str],
                              csv_data: str,
                              distrib: dict,
                              ancestor: int,
                              leaves_only: bool = True,
                              auth: str = "Guest",
                              ) -> dict:
    """Marginalises on an ancestral node using the latent 
    distributions determined from LearnLatentDistributions().

    ONLY implemented for continuous values not discrete multinomial 
    symbols.

    Although its possible, I have not added parameters for 
    rate, seed or gamma values. 

    Parameters:
        nwk(str) = path to file name of nwk 
        states(list) = a list of names for states
        csv_data(str) = path to csv with data
        distrib(dict) = a previously trained distribution from data 
        ancestor(int) = Specify which ancestor to marginalise on
        leaves_only(bool) = ...
        auth(str) = Authentication token, defaults to Guest

    Returns:
        str: {"Message":"Queued","Job":<job-number>}
    """

    request = dict()

    request["Command"] = "Infer"
    request["Auth"] = auth

    params = dict()

    params["States"] = states
    params["Inference"] = "Marginal"
    params["Ancestor"] = ancestor
    params["Leaves-only"] = leaves_only
    params["Distrib"] = distrib

    # format tree
    with open(nwk, 'r') as f:
        tree = ""
        for line in f:
            tree += line.strip()

    params["Tree"] = parsers.nwkToJSON(tree)

    params["Dataset"] = parsers.csvDataToJSON(csv_data)

    request["Params"] = params

    return send_and_recieve(request)


def TrainModes(
        nwk: str,
        csv_data: str,
        n_latent: int,
        leaves_only: bool = None,
        gamma: float = None,
        rounds: int = None,
        seed: float = None,
        rate: float = None,
        alpha: str = "Protein with gap",   # NOTE: currently only works with alphabets pre-defined on Java side
        auth: str = "Guest"
):
    """Provides observations for a given tree and set of evolutionary 'modes'. Operations on the server side optimise
    parameters of modes given input data and tree topology.

    TODO: Much of this is hard-coded for the problem of single-mode problems with motif data and will need to be generalised to work with other data types

    Parameters:
        nwk(str) = Path to nwk file
        csv_data(str) = Path to csv data file
        n_latent(float) = Number of latent states associated with the single mode
        leaves_only(bool) = Consider only plates associated with leaf nodes for training?
        gamma(float) = Uniform substitution model rate parameter
        rounds(int) = Number of rounds for training
        seed(float) = Random seed for training
        rate(float) = Unsure what this rate param is (not JC gamma value?)
        alpha(str) = Domain (alphabet) of characters taken by real nodes
        auth(str) = Authentication token, default is Guest

    """

    request = dict()

    request["Command"] = "TrainModes"
    request["Auth"] = auth

    params = dict()

    params["Dataset"] = parsers.modeDataToJSON(csv_data)

    # format tree
    with open(nwk, 'r') as f:
        tree = ""
        for line in f:
            tree += line.strip()

    params["Tree"] = parsers.nwkToJSON(tree)
    params["Distrib"] = parsers.motifDistribJSON(params["Dataset"]["Features"], n_latent, alpha)

    # If training data contains observations at ancestral nodes, this should be set to False (default True on server side)
    if leaves_only:
        params["Leaves-only"] = leaves_only

    # Optional uniform substitution rate parameter, default set=1.0 on server side
    if gamma:
        params["Gamma"] = gamma

    # Number of rounds for training, default set=10 on server side
    if rounds:
        params["Rounds"] = rounds

    # Random seed for training, default set=1L on server side
    if seed:
        params["Seed"] = seed

    # Unsure what this rate is - appears to be distinct from the JC gamma value. Default is 1 on server side
    if rate:
        params["Rate"] = rate

    request["Params"] = params

    return send_and_recieve(request)


def InferModesMarginal(
        nwk: str,
        csv_data: str,
        distrib: dict,
        query_nodes: list[str],
        infer_latent: bool = None,
        leaves_only: bool = False,
        gamma: float = None,
        rounds: int = None,
        seed : float = None,
        rate : float = None,
        auth: str = "Guest"

):
    """Given distribution parameters for a set of evolutionary modes learnt by TrainModes, infer the marginal
    probability distributions for each feature

    TODO: Much of this is hard-coded for single-mode problems with motif data and will need to be generalised to work with other data types
    E.g. only marginal inference is exposed here (nominated queries are marginally inferred separately)

    Parameters:
        nwk(str) = Path to nwk file
        csv_data(str) = Path to csv data file
        distrib(dict) =
        query_nodes(list) = List of queries for marginal reconstruction (label if leaf, else ancestor number)
        infer_latent(bool) = Should latent nodes be targeted for inference?
        leaves_only(bool) = Consider only plates associated with leaf nodes for inference?
        gamma(float) = Uniform substitution model rate parameter
        rounds(int) = Number of rounds for training
        seed(float) = Random seed for training
        rate(float) = NOTE: Unsure what this rate param is (not JC gamma value?)
        auth(str) = Authentication token, default is Guest

    """

    request = dict()

    request["Command"] = "InferModes"
    request["Auth"] = auth

    params = dict()

    params["Dataset"] = parsers.modeDataToJSON(csv_data)
    params["Distrib"] = distrib

    # format tree
    with open(nwk, 'r') as f:
        tree = ""
        for line in f:
            tree += line.strip()

    params["Tree"] = parsers.nwkToJSON(tree)

    params["Queries"] = query_nodes

    # We'll often be wanting to infer ancestors, so default is set to False here (default = True on server side)
    params["Leaves-only"] = leaves_only

    # TODO: Is this whether to infer the latent state distribution(s) vs. inferring real nodes directly? Unsure
    if infer_latent:
        params["Latent"] = infer_latent

    # Optional uniform substitution rate parameter, default set=1.0 on server side
    if gamma:
        params["Gamma"] = gamma

    # Number of rounds for training, default set=10 on server side
    if rounds:
        params["Rounds"] = rounds

    # Random seed for training, default set=1L on server side
    if seed:
        params["Seed"] = seed

    # Unsure what this rate is - appears to be distinct from the JC gamma value. Default is 1 on server side
    if rate:
        params["Rate"] = rate

    # TODO: Note: Hard-coding for marginal reconstruction here
    params["Inference"] = "Marginal"

    request["Params"] = params

    return send_and_recieve(request)