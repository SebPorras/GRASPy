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
from collections import OrderedDict

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

                tmp["Name"] = key
                tmp["Seq"] = [letter for letter in aln]

                seqs.append(tmp)
            else:
                line = fa.readline()

    alignments = dict()
    alignments["Sequences"] = seqs
    alignments["Datatype"] = {"Predef": "Protein"}

    return alignments


def find_p(s: str):

    start = s.find('(')

    end = s[::-1].find(')')

    if start == -1 and end == -1:
        return start, end

    true_end = len(s) - end - 1

    return start, true_end


def find_comma(s: str, level: int = 0) -> list[str]:

    mylevel = 0
    coms = []

    for i in range(len(s)):

        if s[i] == '(':
            mylevel += 1

        elif s[i] == ')':
            mylevel -= 1
        elif s[i] == ',' and mylevel == level:
            coms.append(i)

    return coms


def nwk_split(n: str, family=None) -> OrderedDict:

    if family is None:
        family = OrderedDict()

    # find the first (children)parent
    start, end = find_p(n)

    # floor case
    if start == -1 and end == -1:

        return family

    # split up the children
    children = n[start+1:end]
    parent = n[end+1:]

    # find every child within the group
    coms = find_comma(children)

    sub_strings = []

    # allows the children to be split up into sub strings
    prev_com = 0

    for i in coms:
        sub_strings.append(children[prev_com:i])
        prev_com = i + 1

    sub_strings.append(children[coms[-1]+1:])

    # this will hold children of the current parent
    current_children = []

    # we must check if the children have children as well
    for i in sub_strings:

        p_start, p_end = find_p(i)

        if p_start == -1 and p_end == -1:

            # if no children, just save them
            current_children.append(i)

        else:
            # break up the children from the paernt
            current_children.append(i[p_end+1:])

            # repeat but with the children of this child
            nwk_split(i, family)

    family[parent] = current_children

    return family


def cleanNwk(nwk: OrderedDict) -> dict:

    # key assumption here is that if there is no distance then it is the root

    dists = dict()
    clean_families = dict()
    labels = []

    for key, values in nwk.items().__reversed__():

        # name:distance
        seq = key.split(':')

        # case for root
        if len(seq) == 1:
            dists[seq[0]] = 0.0
            labels.append(seq[0])

        else:
            dists[seq[0]] = float(seq[1])

        children = []

        for child in values:

            # assign distance
            child = child.split(':')
            dists[child[0]] = float(child[1])

            # record children of current key
            children.append(child[0])

            # add to labels in order
            labels.append(child[0])

        clean_families[seq[0]] = children

    clean_data = {}
    clean_data["Families"] = clean_families
    clean_data["Distances"] = dists
    clean_data["Labels"] = labels

    return clean_data


def nwkToJSON(nwk: str) -> dict:

    #remove ;
    job = nwk[:-1]

    families = nwk_split(job)

    clean = cleanNwk(families)

    idxs = {name: i for i, name in enumerate(clean["Labels"])}

    # PARENTS
    parents = {}

    for key, value in clean["Families"].items():

        for child in value:
            parents[child] = idxs[key]

    for key in clean["Families"].keys():

        if key not in parents.keys():
            parents[key] = -1

    sorted_parents = [int(parents[id]) for id in clean["Labels"]]
    sorted_distances = [clean["Distances"][id] for id in clean["Labels"]]

    json_idx = dict()
    json_idx["Parents"] = sorted_parents
    json_idx["Labels"] = clean["Labels"]
    json_idx["Distances"] = sorted_distances
    json_idx["Branchpoints"] = len(clean["Labels"])

    return json_idx


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

    params["Tree"] = nwkToJSON(tree)
    params["Alignment"] = serialiseAln(aln)

    params["Inference"] = "Joint"
    params["Indels"] = indels
    params["Model"] = model

    request["Params"] = params

    j_request = json.dumps(request) + '\n'

    return client.sendRequest(j_request)


# tree = requestJointReconstruction(aln="./test_data/small_test_data/test_aln.aln",
#                       nwk="./test_data/small_test_data/test_nwk.nwk")

requestPlaceInQueue(27)

output = requestJobOutput(27)
