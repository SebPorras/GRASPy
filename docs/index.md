# GRASPy

GRASPy is the Python API for the Graphical Representation of Ancestral Sequence Predictions.

Read the GRASP paper [here](https://doi.org/10.1371/journal.pcbi.1010633)

## **Requests**

Retrieves any information about a particular job or the output from a job.

### **PlaceInQueue**

    g_requests.PlaceInQueue(job_id:str)

Requests the place in queue of a submitted job.

**Parameters:**

- job_id(str): The ID of the job

**Returns:**

```
str: {"Job": job-number, "Result":{result-JSON}}
```

**Example**

```console

>>> g_requests.requestPlaceInQueue(job_id=19)

Socket created...

Connecting to server...

Socket connected to 10.139.1.21 on IP 4072

Closing socket...

{"Job":19,"Place":0}
```

### **JobOutput**

    g_requests.JobResult(job_id: str)

Requests the output of a submitted job. Request will be
denied if the job is not complete.

**Parameters:**

- job_id(str): The ID of the job

**Returns:**

```
str: {"Job": job-number, "Result":{RESULT}}
```

**Example**

```console

>>>g_requests.JobResult(19)

Socket created...

Connecting to server...

Socket connected to 10.139.1.21 on IP 4072

Closing socket...

```

### **ViewQueue**

    g_requests.ViewQueue()

Lists all the jobs currently being
performed by the server

Parameters:

- job_id(str): The ID of the job

**Returns:**

```
str: {"Job": job-number, "Result":{result-JSON}}
```

**Example**

```console

>>>g_requests.JobResult(19)

Socket created...

Connecting to server...

Socket connected to 10.139.1.21 on IP 4072

{'Jobs': [{'Status': 'COMPLETED','Threads': 1, 'Command': 'Recon','Priority': 0,
'Memory': 1, 'Auth': 'Guest','Job': 1, 'Place': 0},{'Status': 'COMPLETED',
'Threads': 1,'Command': 'Recon','Priority': 0,'Memory': 1, 'Auth': 'Guest',
'Job': 2,'Place': 0}]}

Closing socket...

```

### CancelJob

    g_requests.CancelJob(job_id: str)

Requests the status of a submitted job.

**Parameters:**

- job_id(str): The ID of the job

**Returns:**

```
str: {"Job": job-number}
```

**Example**

```console

>>>g_requests.CancelJob(19)

Socket created...

Connecting to server...

Socket connected to 10.139.1.21 on IP 4072

{'Job': 19}

Closing socket...

```

### **JobStatus**

    g_requests.JobStatus(job_id: str)

Requests the status of job as either completed or queued.

**Parameters:**

- job_id(str): The ID of the job

**Returns:**

```
str: {'Status': 'COMPLETED', 'Job': job-number}
```

**Example**

```console

>>>g_requests.JobStatus(19)

Socket created...

Connecting to server...

Socket connected to 10.139.1.21 on IP 4072

Closing socket...

{'Status': 'COMPLETED', 'Job': job-number}
```

## Commands

### JointReconstruction

    g_requests.JointReconstruction(aln: str, nwk: str, auth: str = "Guest,
    indels: str = "BEP", model: str = "JTT") -> str:

Queries the bnkit server for a joint reconstruction.
Will default to standard bnkit reconstruction parameters which
use BEP for indel model and JTT for the substitution model.

**Parameters:**

- aln(str) = file name of aln file
- nwk(str) = file name of nwk file
- auth(str) = Authentication token, defaults to Guest
- indels(str) = Indel mode, defaults to BEP
- model(str) = Substitution model, defaults to JTT

**Returns:**

```
str - {"Message":"Queued","Job": job-number}
```

**Example**

```console
>>> JointReconstruction(aln="test_aln.aln", nwk="test_nwk.nwk")

Socket created...

Connecting to server...

Socket connected to 10.139.1.21 on IP 4072

Closing socket...

{"Message":"Queued","Job": job-number}

```

### **ExtantPOGTree**

    g_requests.ExtantPOGTree(aln: str, nwk: str, auth: str = "Guest")

Queries the server to turn an alignment
and a nwk file into the POGTree format with POGraphs for extants.
This output JSON can be converted into a POGTree object
via POGTreeFromJSON()

**Parameters:**

- aln(str) = file name of aln file
- nwk(str) = file name of nwk file

**Returns:**

```
Dict: POGTree and POGraphs in JSON format
```

**Example**

```console
>>> ExtantPOGTree(aln="test_aln.aln", nwk="test_nwk.nwk")

Socket created...

Connecting to server...

Socket connected to 10.139.1.21 on IP 4072

Closing socket...

{"Job":<job-number>, "Result":{result-JSON}}

```

### **LearnLatentDistributions**

    LearnLatentDistributions(nwk: str, states: list[str], csv_data: str, auth: str = "Guest")

Learns the distribution of an arbitrary number of discrete
states. The output from the job will be a new/refined distribution.

**CSV input requirements:**

_-CSV files containing data MUST be formatted as outlined below.
The column with names of extant sequences must be named "Headers"_

_-The column with data points MUST be named "Data"._

_- Data with multiple observations must be spaced with whitespace e.g "8.8 1.2"_

```
Headers,Data
A5ILB0,8.5
P08144,7.35 3.3
H9B4I9
```

**Parameters:**

- nwk(str) = path to file name of nwk
- states(list) = a list of names for states
- csv_data(str) = path to csv with data
- auth(str) = Authentication token, defaults to Guest

**Returns:**

```
str: {"Message":"Queued","Job": job-number}
```

**Example**

```console
>>> LearnLatentDistributions(nwk="training.nwk", states=["A", "B"], data="train_data.csv")

Socket created...

Connecting to server...

Socket connected to 10.139.1.21 on IP 4072

Closing socket...

{'Message': 'Queued', 'Job': 42}
```

**Completed job output:**

_- The distributions for each Condition(e.g. "A" or "B") are ordered based on the "Condition" list._

_- Lists within "Pr" contain the mean and variance in that order for each condition._

```console
{ "Distrib":
    { "Condition":[["A"],["B"]],
      "Pr":[[3.784926135903969,0.056738891699391655],
      [2.5324588293595744,0.056738891699391655]],
      "Index":[0,1],
      "Domain":"dat.Continuous@3bd5adde"
    }
}
```

### **MarginaliseDistOnAncestor**

    MarginaliseDistOnAncestor(nwk: str, states: list[str], csv_data: str, distrib: dict, ancestor: int, leaves_only: bool = True, auth: str = "Guest"):

Marginalises on an ancestral node using the latent
distributions determined from LearnLatentDistributions().
Although its possible, I have not added parameters for
rate, seed or gamma values.

**CSV input requirements:**

_-CSV files containing data MUST be formatted as outlined below.
The column with names of extant sequences must be named "Headers"_

_-The column with data points MUST be named "Data"._

_- Data with multiple observations must be spaced with whitespace e.g "8.8 1.2"_

```
Headers,Data
A5ILB0,8.5
P08144,7.35 3.3
H9B4I9
```

**Parameters:**

- nwk(str) = path to file name of nwk
- states(list) = a list of names for states
- csv_data(str) = path to csv with data
- distrib(dict) = a previously trained distribution from data
- ancestor(int) = Specify which ancestor to marginalise on
- leaves_only(bool) = ...
- auth(str) = Authentication token, defaults to Guest

**Returns:**

```
str: {"Message":"Queued","Job": job-number}
```

**Example:**

```console
>>> MarginaliseDistOnAncestor(nwk="training.nwk", states=["A", "B"], data="train_data.csv", ancetor=0)

Socket created...

Connecting to server...

Socket connected to 10.139.1.21 on IP 4072

Closing socket...

{'Message': 'Queued', 'Job': 42}
```

**Completed job output:**

_- The marginalised distribution for each state are ordered in the list according to "Values"._

_- Values within "Pr" represent the mean and variance respectively._

```console
{ "N0":[
    { "Pr":[0.6652270145537978,0.3347729854462022],
      "Domain":{"Size":2,"Values":["A","B"],"Datatype":"String"}},
    { "Pr":[0.649968113685095,0.350031886314905],
    "Domain":{"Size":2,"Values":["A","B"],"Datatype":"String"}}
    ]
}
```

## **Data Structures**

- GRASPy has a number of data structures that can be used to interact with output from the bnkit server responses.

### **POGTreeFromJointReconstruction**

    parsers.POGTreeFromJSON(nwk: Union[str, dict], POG_graphs: dict)

Creates an instance of the POGTree data structure. A nwk
file OR output from JointReconstruction() can be used
to create tree topology with the second option also creating
POGraphs for extants.

**Parameters:**

- nwk (str or dict): Users can input a nwk file path or can provide
  the output from g_requests.requestPOGTree().
- POG_graphs(dict): The POGraphs for ancestors generated from
  output from g_requests.requestJointReconstruction().

**Returns:**

```
<POGTree>
```

**Example:**

```console
>>> tree = POGTreeFromJointReconstruction(nwk="example.nwk", POG_graphs=graphs)
```

### **POGTree**

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

**Parameters:**

- idxTree: Instance of the IdxTree class
- POGraphs(dict): maps sequence IDs to POGraph class
- nBranches(int): number of branch points in the tree
- branchpoints(dict[str, BranchPoint]): Contains BranchPoint objects
- parents(list[int]): maps the index of the child to the index
  of the parent
- children(list[list[Union[int, None]]]): maps the index of the
  parent to an array containing the indexes the child/children
- indices(dict[str, int]): Maps the sequence ID to the index on the
  tree
- distances(list[float]): Maps the branchpoint to the distance to
  its parent

#### **writeToNwk**

    writeToNwk(file_name: str, root: str = "N0")

Converts the POGTree object into a nwk string and writes this to a file

**Parameters:**

- file_name(str) : name of nwk file
- root(str): Default set to N0 at the "root" ancestor but can be changed to internal nodes to create subtrees if desired.

**Returns:**

```
str: The POGTree in nwk format
```

**Example:**

```console
>>> tree = POGTree(nwk.nwk, aln.aln)
>>> nwk = tree.writeNwk(test_nwk.nwk)
>>> print(nwk)
(XP_004050792.2:0.040380067,XP_005216113.1:0.028035396,(XP_018963554.1:0.016721581,XP_016357833.1:0.024301326)N1:0.347992941)N0:0;
```

#### **writeToFasta**

    writeToFasta(file_name: str)

Writes all sequences of the tree to file.
Sequence for ancestors are based on a joint
reconstruction and each symbol is the most likely
at each position.

**Parameters:**

- file_name(str) : name of fasta file

**Returns:**

```
None
```

**Example:**

```console
>>>POGTree("test_fasta")
>>>
```

### **BranchPoint**

    BranchPoint(id: str, parent: Union[str, None], dist: float,
                 children: list[str], seq: Optional[Sequence] = None)

Represents a branchpoint on a phylogenetic tree. Can contain
information about the parents or children of that point and how
long that branch point is.

**Parameters:**

- id(str): Sequence ID
- parent(str or None): ID of parent
- dist(float): Distance to parent
- children(list): IDs children of current BranchPoint
- seq(Sequence): The sequence based on a joint reconstruction
  if the BranchPoint is an ancestor otherwise it is just
  the sequence of an extant. Contains a Sequence object.

### **SymNode**

    SymNode(name: int, symbol: str, edges: list)

Only implemented for output from joint reconstruction.
Stores the most likely character at a given sequence position
and all of the edges at this position.

**Parameters:**

- name(int): index position in sequence
- symbol(str): Most likely amino acid based on joint reconstruction
- edges(list): Contains all outgoing edges at this position

### **Edge**

    Edge(start: int, end: int, edgeType: Optional[str] = None, recip: Optional[bool] = None, backward: Optional[bool] = None, forward: Optional[bool] = None, weight: Optional[float] = None)

Creates instance of an edge between two positions in a sequence.
Currently only implemented for bidirectional edges.

**Parameters:**

- start(int): position of beginning of edge
- end(int): position of end of edge
- edgeType(str): Currently only supports bidirectional edge
- recip(bool): ASK ABOUT THIS
- backward(bool): Direction of edge
- forward(bool): Direction of edge
- weight(float): Support of the edge
