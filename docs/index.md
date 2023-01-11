# GRASPy 

GRASPy is the python API for the Graphical Representation of Ancestral Sequence 
Predictions

Read the paper [here](https://doi.org/10.1371/journal.pcbi.1010633)

## Commands 
- GRASPy interacts with our servers via sockets. The following commands can be 
used to submit and retrieve your jobs.

### requestJointReconstruction

    commands.requestJointReconstruction(aln: str, nwk: str, auth: str = "Guest,
    indels: str = "BEP", model: str = "JTT") -> str:

Queries the bnkit server for a joint reconstruction.
Will default to standard bnkit reconstruction parameters which
use BEP for indel model and JTT for the substitution model.

Parameters:

> * aln(str) = file name of aln file
> * nwk(str) = file name of nwk file
> * auth(str) = Authentication token, defaults to Guest
> * indels(str) = Indel mode, defaults to BEP
> * model(str) = Substitution model, defaults to JTT

Returns: 
> str - returns Job ID

Example

```console 
>>> requestJointReconstruction(aln="test_aln.aln", nwk="test_nwk.nwk")

Establishing connectoin...

Request sent...

Your job ID is: 123ABC
```

### requestJobStatus

    commands.requestJobStatus(job_id: str, auth: str = "Guest")

Requests the status of a submitted job. 

Parameters:

> * job\_id(str): The ID of the job 
> * auth(str): Authentication token, defaults to Guest

Returns: 

> str: returns job status 

Example

```console 
>>> requestJobStatus(job_id = "ABC123", auth = "Guest")

Establishing connection...

Request sent...

Your job status is: <Job_status>
```

### requestJobResult 

    commands.requestJobResult(job_id: str, auth: str = "Guest")

Requests the output of a submitted job. Request will be 
denied if the job is not complete. 

Parameters:

> * job\_id(str): The ID of the job 
> * auth(str): Authentication token, defaults to Guest

Returns:

> str: returns completed job 

Example

```console 
>>> requestJobResult(job_id = "ABC123", auth = "Guest")

Establishing connection...

Request sent...

Your job result is: <job_output>
```

## Data Structures

* GRASPy has a number of data structures that can be used to interact with output from the bnkit server responses. 


### POGTreeFromJSON

    POGTreeFromJSON(json_path: str)

Static function that can instantiate a POGTree object from a JSON file.

Parameters:

> * json\_path(str): path to JSON file  

Returns: 
> POGTree

Example

```console 
>>> tree = POGTreeFromJSON("ASR.json")
```

### POGTree
    
    POGTree(idxTree: IdxTree, POGraphs: dict[str, POGraph])

The Partial Order Graph Tree (POGTree), is a phylogenetic tree that contains two data structures. 1) An IdxTree object which holds information about the topology and information about the sequence at each branchpoint of the tree. 2) A POGraph object which describes the graph of the sequence at that branchpoint.

Parameters:

> * idxTree: Instance of the IdxTree class  
> * POGraphs(dict): maps sequence IDs to POGraph class  


#### writeNwk 

    pog_tree.writeNwk(file_name: str, root: str = "N0")

Converts the POGTree object into a nwk string and writes this to a file 

Parameters:

> * file\_name(str) : name of nwk file 
> * root(str): Default set to N0 at the "root" ancestor but can be changed to internal nodes to create subtrees if desired.

Returns: 
> str: The nwk string that was written to the file

Example

```console 
>>> tree = POGTree(nwk.nwk, aln.aln)
>>> nwk = tree.writeNwk(test_nwk.nwk)
>>> print(nwk)
(XP_004050792.2:0.040380067,XP_005216113.1:0.028035396,(XP_018963554.1:0.016721581,XP_016357833.1:0.024301326)N1:0.347992941)N0:0;
```


#### template  

    pog_tree.writeNwk(file_name: str, root: str = "N0")

Converts the POGTree object into a nwk string and writes this to a file 

    Parameters:

file\_name(str) : name of nwk file 
    
root(str): Default set to N0 at the "root" ancestor but can be changed to internal nodes to create subtrees if desired.

Returns: 
    str: The nwk string that was written to the file

Example

```console 
>>> tree = POGTree(nwk.nwk, aln.aln)
>>> nwk = tree.writeNwk(test_nwk.nwk)
>>> print(nwk)
(XP_004050792.2:0.040380067,XP_005216113.1:0.028035396,(XP_018963554.1:0.016721581,XP_016357833.1:0.024301326)N1:0.347992941)N0:0;
```

### IdxTree

    IdxTree( nBranches: int, branchpoints: dict[str, BranchPoint], parents: list[int], children: list[list[Union[int, None]]], indices: dict[str, int], distances: list[float])

IdxTree is condensed phylogenetic tree representation within a POGTree object. Each branchpoint is assigned an index and a BranchPoint object, allowing easy access of information via the sequence name of an extant or an ancestor. 
