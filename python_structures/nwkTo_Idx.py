import json

nwk = "(((XP_006629927.2:0.14777721590070114,(XP_018611667.1:0.31741219299255685,((XP_007229530.1:0.20410994427385432,(XP_012687241.1:0.12553661262814741,XP_018919739.1:0.17752015624894124)N6:0.02108546966591862)N5:0.024498458258308764,(XP_019717376.1:0.3374131974537633,(XP_010886716.1:0.17046994270732152,(XP_014050304.1:0.012222797694125376,XP_021429054.1:0.016443894223391986)N9:0.051214958442506786)N8:0.06898824252200719)N7:0.07186524457989396)N4:0.10310490468352707)N3:0.06898460428141462)N2:0.11632291055872579,(ARO89866.1:0.25205280271933694,(XP_014733783.1:0.29681349725859363,(XP_005082857.1:0.1589210754372825,((XP_021540185.1:0.03325108724429082,XP_019684690.2:0.018088332404641694)N14:0.04823406553019316,(XP_012621711.1:0.05681896074666115,((NP_898898.1:1.0000005001842283e-06,XP_526649.2:0.00220424589188406)N17:0.005288422325157249,(XP_012291909.1:0.0046099502113805535,XP_003929520.1:0.004241068270090409)N18:0.010508512073819531)N16:0.040465306764365216)N15:0.031854968122094096)N13:0.030935504331779384)N12:0.18757970335837193)N11:0.08013548717836616)N10:0.11264351323303923)N1:0.7819801484223995,((XP_004050792.2:0.046745089073882085,XP_005216113.1:0.031339011068860056)N20:0.22157549651088937,(XP_018963554.1:0.019918423012584174,XP_016357833.1:0.026700451827796012)N21:0.2239611736832119)N19:0.7819801484223997)N0:0;"


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


def nwk_split(n: str, family=None) -> dict:

    if family is None:
        family = {}

    start, end = find_p(n)

    if start == -1 and end == -1:

        return family

    children = n[start+1:end]
    parent = n[end+1:]

    coms = find_comma(children)

    sub_strings = []

    prev_com = 0

    for i in coms:
        sub_strings.append(children[prev_com:i])
        prev_com = i + 1

    sub_strings.append(children[coms[-1]+1:])

    sub_family = []

    for i in sub_strings:

        p_start, p_end = find_p(i)

        if p_start == -1 and p_end == -1:

            sub_family.append(i)

        else:
            # adds the head of the new group
            sub_family.append(i[p_end+1:])

            # repeat but with the subgroup
            nwk_split(i, family)

    family[parent] = sub_family

    return family


def nwkToJSON(nwk: str) -> str:

    families = nwk_split(nwk)

    # LABELS AND DISTS

    labels = []
    distances = []

    # # reverse to go in order from root to leaves
    for key in families.keys().__reversed__():

        labels.append(key.split(':')[0])
        distances.append(float(key.split(':')[1].replace(';', '')))

    idxs = {name: i for i, name in enumerate(labels)}

    # PARENTS
    parents = {}

    for key, value in families.items():

        for child in value:
            parents[child.split(':')[0]] = idxs[key.split(':')[0]]

    for key in families.keys():

        if key.split(':')[0] not in parents.keys():
            parents[key.split(':')[0]] = -1

    sorted_parents = [int(parents[id]) for id in labels]

    json_idx = dict()
    json_idx["Parents"] = sorted_parents
    json_idx["Labels"] = labels
    json_idx["Distances"] = distances
    json_idx["Branchpoints"] = len(labels)

    j_format = json.dumps(json_idx)

    return j_format


test = nwkToJSON(nwk)
