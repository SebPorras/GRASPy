# GRASPy

GRASPy is the python API for the Graphical Representation of Ancestral Sequence
Predictions

Read the paper [here](https://doi.org/10.1371/journal.pcbi.1010633)

## Requests

- GRASPy interacts with our servers via sockets. The following commands can be
  used to submit and retrieve your jobs.

### requestPlaceInQueue

    g_requests.requestPlaceInQueue(job_id:str)

Requests the place in queue of a submitted job.

Parameters:

> - job_id(str): The ID of the job

Returns:

> str: {"Job":<job-number>, "Result":{<result-JSON>}}

Example

```console

>>> g_requests.requestPlaceInQueue(job_id=19)

Socket created...

Connecting to server...

Socket connected to 10.139.1.21 on IP 4072

Closing socket...

'{"Job":19,"Place":0}
```

### requestJobOutput

    g_requests.requestJobResult(job_id: str)

Requests the output of a submitted job. Request will be
denied if the job is not complete.

Parameters:

> - job_id(str): The ID of the job

Returns:

> str: {"Job":<job-number>, "Result":{RESULT}}

Example

```console

>>>g_requests.requestJobResult(19)

Socket created...

Connecting to server...

Socket connected to 10.139.1.21 on IP 4072

Closing socket...

```

### requestViewQueue

    g_requests.requestViewQueue()

Lists all the jobs currently being
performed by the server

Parameters:

> - job_id(str): The ID of the job

Returns:

> str: {"Job":<job-number>, "Result":{<result-JSON>}}

Example

```console

>>>g_requests.requestJobResult(19)

Socket created...

Connecting to server...

Socket connected to 10.139.1.21 on IP 4072

{'Jobs': [{'Status': 'COMPLETED','Threads': 1, 'Command': 'Recon','Priority': 0,
'Memory': 1, 'Auth': 'Guest','Job': 1, 'Place': 0},{'Status': 'COMPLETED',
'Threads': 1,'Command': 'Recon','Priority': 0,'Memory': 1, 'Auth': 'Guest',
'Job': 2,'Place': 0}]}

Closing socket...

```

### requestCancelJob

    g_requests.requestCancelJob(job_id: str)

Requests the status of a submitted job.

Parameters:

> - job_id(str): The ID of the job

Returns:

> str: {"Job":<job-number>}

Example

```console

>>>g_requests.requestCancelJob(19)

Socket created...

Connecting to server...

Socket connected to 10.139.1.21 on IP 4072

{'Job': 19}

Closing socket...

```

## Commands

### requestJointReconstruction

    g_requests.requestJointReconstruction(aln: str, nwk: str, auth: str = "Guest,
    indels: str = "BEP", model: str = "JTT") -> str:

Queries the bnkit server for a joint reconstruction.
Will default to standard bnkit reconstruction parameters which
use BEP for indel model and JTT for the substitution model.

Parameters:

> - aln(str) = file name of aln file
> - nwk(str) = file name of nwk file
> - auth(str) = Authentication token, defaults to Guest
> - indels(str) = Indel mode, defaults to BEP
> - model(str) = Substitution model, defaults to JTT

Returns:

> str - Job number

Example

```console
>>> requestJointReconstruction(aln="test_aln.aln", nwk="test_nwk.nwk")

Socket created...

Connecting to server...

Socket connected to 10.139.1.21 on IP 4072

Closing socket...

{"Message":"Queued","Job":19}

```

### requestExtantPOGTree

    g_requests.requestExtantPOGTree(aln: str, nwk: str, auth: str = "Guest")

Queries the server to turn an alignment
and a nwk file into the POGTree format with POGraphs for extants.
This output JSON can be converted into a POGTree object
via POGTreeFromJSON()

Parameters:

> - aln(str) = file name of aln file
> - nwk(str) = file name of nwk file

Returns:

> Dict: POGTree and POGraphs in JSON format

Example

```console
>>> requestJointReconstruction(aln="test_aln.aln", nwk="test_nwk.nwk")

Socket created...

Connecting to server...

Socket connected to 10.139.1.21 on IP 4072

Closing socket...

{"Job":<job-number>, "Result":{<result-JSON>}}

```

### requestTrainFromData

    requestTrainFromData(nwk: str, states: list[str], data: str, auth: str "Guest"):

Learns the distribution of an arbitrary number of discrete
states. The output from the job will be a new/refined distribution.

Parameters:

> - nwk(str) = path to file name of nwk
> - states(list) = a list of names for states
> - data(str) = path to csv with data
> - auth(str) = Authentication token, defaults to Guest

Returns:

> str - Job number

Example

```console
>>> requestTrainFromData(nwk="./3_2_1_1_filt.nwk", states=["A", "B"], data="train_data.csv")

Socket created...

Connecting to server...

Socket connected to 10.139.1.21 on IP 4072

Closing socket...

{"Job":<job-number>, "Result":{<result-JSON>}}
```

## Data Structures

- GRASPy has a number of data structures that can be used to interact with output from the bnkit server responses.

### POGTreeFromJSON

    parsers.POGTreeFromJSON(nwk: Union[str, dict], POG_graphs: dict)

Creates an instance of the POGTree data structure. A nwk
file OR output from g_requests.requestPOGTree() can be used
to create tree topology with the second option also creating
POGraphs for extants.

Parameters:

> - nwk (str or dict): Users can input a nwk file path or can provide
>   the output from g_requests.requestPOGTree().
> - POG_graphs(dict): The POGraphs for ancestors generated from
>   output from g_requests.requestJointReconstruction().

Returns:

> POGTree

Example

```console
>>> tree = parsers.POGTreeFromJSON(nwk="example.nwk", POG_graphs=graphs)
```

### POGTree

    POGTree(nBranches: int, branchpoints: dict[str, BranchPoint],
                 parents: list[int], children: list[list[Union[int, None]]],
                 indices: dict[str, int], distances: list[float],
                 POGraphs: dict[str, POGraph])

The Partial Order Graph Tree (POGTree), is a phylogenetic tree made up
of branchpoints which represent nodes on the tree.
Each branchpoint is assigned an index and a BranchPoint
object, allowing easy access of information via the sequence name of
an extant or an ancestor. 2) A POGraph object which describes
the graph of the sequence at that branchpoint.

Parameters:

> - idxTree: Instance of the IdxTree class
> - POGraphs(dict): maps sequence IDs to POGraph class
> - nBranches(int): number of branch points in the tree
>   branchpoints(dict[str, BranchPoint]): Contains BranchPoint objects
> - parents(list[int]): maps the index of the child to the index
>   of the parent
> - children(list[list[Union[int, None]]]): maps the index of the
>   parent to an array containing the indexes the child/children
> - indices(dict[str, int]): Maps the sequence ID to the index on the
>   tree
> - distances(list[float]): Maps the branchpoint to the distance to
>   its parent

POGraph(POGraph): Instance of a POGraph object

#### writeToNwk

    pog_tree.writeToNwk(file_name: str, root: str = "N0")

Converts the POGTree object into a nwk string and writes this to a file

Parameters:

> - file_name(str) : name of nwk file
> - root(str): Default set to N0 at the "root" ancestor but can be changed to internal nodes to create subtrees if desired.

Returns:

> str: The POGTree in nwk format

Example

```console
>>> tree = POGTree(nwk.nwk, aln.aln)
>>> nwk = tree.writeNwk(test_nwk.nwk)
>>> print(nwk)
(XP_004050792.2:0.040380067,XP_005216113.1:0.028035396,(XP_018963554.1:0.016721581,XP_016357833.1:0.024301326)N1:0.347992941)N0:0;
```
