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


def serialiseAln(file_name: str) -> list[dict]:
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

    return seqs


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
        str: returns Job ID
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
    params["Alignments"] = serialiseAln(aln)
    params["Inference"] = "Joint"
    params["Indels"] = indels
    params["Model"] = model

    request["Parameters"] = params

    j_request = json.dumps(request) + '\n'
    print("j_request")

    print(j_request)

    return client.sendRequest(j_request)


def requestJobStatus(job_id: str, auth: str = "Guest") -> str:
    """Requests the status of a submitted job

    Parameters:
        job_id(str): The ID of the job 
        auth(str): Authentication token, defaults to Guest

    Returns:
        str: bnkit job status  
    """

    request = dict()

    request["Command"] = "JobStatus"
    request["Auth"] = auth
    request["Job"] = job_id

    j_request = json.dumps(request) + '\n'

    return client.sendRequest(j_request)


def requestJobResult(job_id: str, auth: str = "Guest") -> str:
    """Requests the output of a submitted job. Request will be 
    denied if the job is not complete. 

    Parameters:
        job_id(str): The ID of the job 
        auth(str): Authentication token, defaults to Guest

    Returns:
        str: Bnkit response for completed job 
    """

    request = dict()

    request["Command"] = "JobResult"
    request["Auth"] = auth
    request["Job"] = job_id

    j_request = json.dumps(request) + '\n'

    return client.sendRequest(j_request)


joint = requestJointReconstruction(
    aln="./test_data/small_test_data/test_aln.aln",
    nwk="./test_data/small_test_data/test_nwk.nwk")

#jobStatus = requestJobStatus("ABC123")

# JobResult = requestJobResult("ABC123")
