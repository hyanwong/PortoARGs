import collections
import itertools
import sys
import numpy as np
from IPython.display import HTML


### Hack below - can be removed when questions are saved to JSON files named Q1.json, Q2.json etc
# This allows answers to be hidden from the casual viewer. Then we can simply set the url to a string instead of a class
class FakeURL(dict):
    def __add__(self, prefix):
        return self[prefix]
WB1_base = FakeURL()
true, false= True, False  # just to allow easy conversion to JSON format
### hack ends

WB1_base["Q1.json"] = [{
    "question": "What genomic span is passed on to sample node 6 through the route via node 15? The span is the total length of the genome region passed on (i.e. the right minus the left position)",
    "type": "numeric",
    "answers": [
        {"type": "value", "value": 930, "correct": true, "feedback": "Yes, it covers positions 0 to 930."},
        {"type": "default", "feedback": "Try hovering over the edges below node 15."}
    ]
}]
WB1_base["Q2.json"] = [{
    "question": "What does SPR stand for?",
    "type": "many_choice",
    "answers": [
        {"answer": "Subtree Prune and Regraft", "correct": true},
        {"answer": "Strand Pairing Resolution", "correct": false},
        {"answer": "Single-Point Reversion", "correct": false},
        {"answer": "Slippery Puzzle Reorganization", "correct": false, "feedback": "???."}
    ]}, {
    "question": "What is the TMRCA, or time to the most recent common ancestor (in generations) between nodes 0 and 9 at the left hand end of this ARG?",
    "type": "numeric",
    "answers": [
        {"type": "value", "value": 38472, "correct": true},
        {"type": "value", "value": 31, "correct": false, "feedback": "That's the node ID. What time is associated with that node."},
        {"type": "default", "feedback": "Try hovering over the last part of the genome, and looking at the axis labels."}
    ]},{
    "question": "What is the name of the structure in which two lineages split but immediately re-join, as seen just below node 26",
    "type": "many_choice",
    "answers": [
        {"answer": "Diamond", "correct": true},
        {"answer": "Loop", "correct": false},
        {"answer": "Bubble", "correct": false, "feedback": "This term is sometimes used, but is not standard"},
        {"answer": "Cycle", "correct": false}
    ]
}]

WB1_base["Q3.json"] = [{
    "question": "(Hard!) Only one recombination event results in a tree that shows a new topological relationship between the samples. By looking at the graph and the local trees, can you identify which one it is? Enter its breakpoint position below:",
    "type": "numeric",
    "answers": [
        {"type": "value", "value": 930, "correct": true, "feedback": "Yes (coincidentally this is the last breakpoint)"},
        {"type": "default", "feedback": "Hint: look at the breakpoints on the X axis of the tree-by-tree plot: which marks the transition between two differently shaped trees?"}
    ]
}]

WB1_base["Q4.json"] = [{
    "question": "The first two breakpoints along the genome happen to correspond to recombination at the top of the ARG. How many different root heights does this cause?",
    "type": "numeric",
    "answers": [
        {"type": "value", "value": 3, "correct": true},
        {"type": "default", "feedback": "Hint: look at the local trees"}
    ]
}]

WB1_base["Q5.json"] = [{
    "question": "How would you classify the recombination node(s) 20/21",
    "type": "many_choice",
    "answers": [
        {"answer": "Non-coalescent / always unary", "correct": true, "feedback": "By construction, in a full ARG recombination nodes are never coalescent"},
        {"answer": "Partially-coalescent / locally unary", "correct": false},
        {"answer": "All-coalescent / never unary", "correct": false}
    ]},{
    "question": "How would you classify node 17",
    "type": "many_choice",
    "answers": [
        {"answer": "Non-coalescent / always unary", "correct": false},
        {"answer": "Part-coalescent / locally unary", "correct": false},
        {"answer": "All-coalescent / never unary", "correct": true}
    ]},{
    "question": "How would you classify node 26",
    "type": "many_choice",
    "answers": [
        {"answer": "Non-coalescent / always unary", "correct": true, "feedback": "Some 'common ancestor' nodes with 2 children in the full ARG never represent local coalescence"},
        {"answer": "Part-coalescent / locally unary", "correct": false},
        {"answer": "All-coalescent / never unary", "correct": false}
    ]},{
    "question": "How would you classify node 15",
    "type": "many_choice",
    "answers": [
        {"answer": "Non-coalescent / always unary", "correct": false},
        {"answer": "Part-coalescent / locally unary", "correct": true, "feedback": "It is only coalescent on the left side of the genome: it is unary to the right of position 930"},
        {"answer": "All-coalescent / never unary", "correct": false}
    ]
}]


WB1_base["Q6.json"] = [{
    "question": "In the population-coloured ARG plot, what colour do you think represents nodes from the African population",
    "type": "many_choice",
    "answers": [
        {"answer": "green", "correct": true, "feedback": "Yes, the deepest divergences are African: in fact, even between the European genomes, some coalescences trace back into Africa"},
        {"answer": "blue", "correct": false},
    ]
}]

WB1_base["Q7.json"] = [{
    "question": "One of the sites has two mutations. Can you identify the position of that site by looking at both the tree-by-tree and the interactive ARG plot?",
    "type": "numeric",
    "answers": [
        {"type": "value", "value": 5, "correct": true, "feedback": "Yes (coincidentally this is the first position)"},
        {"type": "default", "feedback": "Hint: duplicate mutations will simultaneously be highlighted when hovering over them in the visualizer"}
    ]
}]

WB1_base["Q8.json"] = [{
    "question": "What is the allelic state of the sample with node ID 2 at position 5?",
    "type": "many_choice",
    "answers": [
        {"answer": "A", "correct": false},
        {"answer": "C", "correct": false},
        {"answer": "G", "correct": false},
        {"answer": "T", "correct": true}
    ]}, {
    "question": "What is the ancestral state at position 5?",
    "type": "many_choice",
    "answers": [
        {"answer": "A", "correct": true},
        {"answer": "C", "correct": false},
        {"answer": "G", "correct": false},
        {"answer": "T", "correct": false}
    ]}, {
    "question": "By looking at the ARG or tree visualizations, the mutation responsible for this state in sample 2 is above which node?",
    "type": "many_choice",
    "answers": [
        {"answer": "12", "correct": false},
        {"answer": "13", "correct": true},
        {"answer": "31", "correct": false},
        {"answer": "This allelic state is the ancestral state, so does not correspond to a mutation in the ARG ", "correct": false,
        "feedback": "It is true that ancestral states are not represented by a mutation in the ARG, but the state C here is a derived state, so it *is* associated with a mutation"}
    ]
}]

WB1_base["Q9.json"] = [{
    "question": "In the ARG above, what mutation groups together samples 7, 8, and 9?",
    "type": "many_choice",
    "answers": [
        {"answer": "C mutates to T at position 215", "correct": true},
        {"answer": "T mutates to C at position 215", "correct": false, "feedback": "You have the position right, but the derived state is T."},
        {"answer": "G mutates to T at position 92", "correct": false},
        {"answer": "C mutates to G at position 560", "correct": false}
    ]}, {
    "question": "What ARG node does that mutation allow us to infer?",
    "type": "numeric",
    "answers": [
        {"type": "value", "value": 16, "correct": true, "feedback": "Yes, it's the node below the mutation."},
        {"type": "value", "value": 19, "correct": false, "feedback": "No, that node groups together more than 7, 8, and 9."}
    ]},{
    "question": "Do any mutations allow us to resolve the relationship between samples 4, 5, and 6?",
    "type": "many_choice",
    "answers": [
        {"answer": "No", "correct": true, "feedback": "In fact, we can't even resolve the relative closeness of 4 versus 5 versus 6 to the (7,8,9) group"},
        {"answer": "Yes", "correct": false}
    ]}, {
    "question": "How many mutations are uninformative about the structure (topology) of the ARG?",
    "type": "numeric",
    "answers": [
        {"type": "value", "value": 2, "correct": true, "feedback": "Yes, the two so-called singletons, immediately above node 3 and node 4."},
        {"type": "default", "feedback": "Incorrect (hint: mutations above a single node do not group samples together)."}
    ]
}]

WB1_base["Q10.json"] = [{
    "question": "What is the interpretation of the genetic diversity, π in terms of branch lengths?",
    "type": "many_choice",
    "answers": [
        {"answer": "The average time of all internal nodes", "correct": false},
        {"answer": "Twice the average TMRCA between all pairs of samples", "correct": true},
        {"answer": "The average branch lengths between all pairs of samples", "correct": true},
        {"answer": "The span-weighted average branch length of the local trees", "correct": false,
        "feedback": "This is the equivalent of the number of segregating sites (under the infinite-sites model)"}
    ]
}]

WB1_base["Q11.json"] = [{
    "question": "In the diamond in our example ARG, what is the ID of the node representing its common ancestor (CA) event?",
    "type": "numeric",
    "answers": [
        {"type": "value", "value": 26, "correct": true},
        {"type": "default"}
    ]
}]

WB1_base["Q12.json"] = [{
    "question": "What colour are all-coalescent nodes in the graph above?",
    "type": "many_choice",
    "answers": [
        {"answer": "Green", "correct": true},
        {"answer": "Blue", "correct": false},
        {"answer": "Red", "correct": false},
        {"answer": "Cyan", "correct": false},
    ]},{
    "question": "What colour are partly-coalescent nodes in the graph above?",
    "type": "many_choice",
    "answers": [
        {"answer": "Green", "correct": false},
        {"answer": "Blue", "correct": true},
        {"answer": "Red", "correct": false},
        {"answer": "Cyan", "correct": false},
    ]},{
    "question": "What colour are non-coalescent CA nodes in the graph above?",
    "type": "many_choice",
    "answers": [
        {"answer": "Green", "correct": false},
        {"answer": "Blue", "correct": false},
        {"answer": "Red", "correct": false, "feedback": "Not quite: red nodes are indeed non-coalescent, but they are RE nodes"},
        {"answer": "Cyan", "correct": true, "feedback": "Correct (the other non-coalescent nodes are the red RE nodes)"},
    ]},{
    "question": "Are non-coalescent CA event nodes always above a diamond?",
    "type": "many_choice",
    "answers": [
        {"answer": "No", "correct": true, "feedback": "Node 35 is not above a diamond, and it is non-coalescent in every local tree in which it is seen."},
        {"answer": "Yes", "correct": false},
    ]
}]

WB1_base["Q13.json"] = [{
    "question": "What has happened to the diamond and the other non-coalescent CA node?",
    "type": "many_choice",
    "answers": [
        {"answer": "They have been bypassed entirely", "correct": true},
        {"answer": "They have been kept", "correct": false},
        {"answer": "They are now connected to other nodes", "correct": false}
    ]}, {
    "question": "Other than those cases, what has happened to the edges descending from the parents of the red recombination nodes?",
    "type": "many_choice",
    "answers": [
        {"answer": "They now bypass the recombination node and link to the next (blue or green) node below", "correct": true},
        {"answer": "They have been deleted entirely, and not replaced", "correct": false},
        {"answer": "They have been shortened", "correct": false}
    ]}, {
    "question": "What is the maximum number of parents of any node in the new ARG?",
    "type": "numeric",
    "answers": [
        {"type": "value", "value": 3, "correct": true, "feedback": "Correct: node 23 now represents where two recombination events have been collapsed into a single node below"},
        {"type": "default"}
    ]}, {
    "question": "What is the maximum number of children of any node in the new ARG?",
    "type": "numeric",
    "answers": [
        {"type": "value", "value": 3, "correct": true, "feedback": "Correct: the root node 36 now represents where two coalescent events have been collapsed into a single node above"},
        {"type": "default"}
    ]
}]

WB1_base["Q14.json"] = [{
    "question": "After completely deleting unknowable nodes, what is the ID of the node with 3 children in the ARG?",
    "type": "numeric",
    "answers": [
        {"type": "value", "value": 22, "correct": true, "feedback": "Correct. However, note that it still only has 2 children in any one local tree"},
        {"type": "value", "value": 36, "correct": false, "feedback": "That's the old node ID: you need to simplify without the filter_nodes=False option"},
        {"type": "default"}
    ]}, {
    "question": "After completely deleting unknowable nodes, what is the ID of the node with 3 parents in the ARG?",
    "type": "numeric",
    "answers": [
        {"type": "value", "value": 19, "correct": true, "feedback": "Correct. However, note that it's still the case that all nodes only have one parent in any local tree"},
        {"type": "value", "value": 22, "correct": false, "feedback": "That's the old node ID: you need to simplify without the filter_nodes=False option"},
        {"type": "default"}
    ]
}]

WB1_base["Q15.json"] = [{
    "question": "In the simplified ARG, do we know the order of the two recombination events immediately above node 19",
    "type": "many_choice",
    "answers": [
        {"answer": "No", "correct": true},
        {"answer": "Yes", "correct": false},
    ]}, {
    "question": "In the simplified ARG, do we know the order of the two coalescent events immediately below the root node 22",
    "type": "many_choice",
    "answers": [
        {"answer": "No", "correct": true},
        {"answer": "Yes", "correct": false},
    ]}]

WB1_base["Q16.json"] = [{
    "question": "Previously, there were 7 local trees corresponding to 6 recombination breakpoints. How many breakpoints are now represented?",
    "type": "numeric",
    "answers": [
        {"type": "value", "value": 5, "correct": true, "feedback": "Correct: the recombination at position 602, associated with the diamond, has been removed."},
        {"type": "default"}
    ]}, {
    "question": "Previously we had 13 mutations at twelve sites. How many mutations are there in the simplified ARG?",
    "type": "numeric",
    "answers": [
        {"type": "value", "value": 13, "correct": true, "feedback": "Correct: simplification doesn't change the sites or the number / age of mutations"},
        {"type": "default"}
    ]
}]

WB1_base["Q17.json"] = [{
    "question": "Compared to the partially simplified ARG, has the number of EDGES in the fully simplified ARG increased, decreased, or stayed the same?",
    "type": "many_choice",
    "answers": [
        {"answer": "Increased", "correct": true, "feedback": "Correct: full simplification can actually increase the number of ARG edges"},
        {"answer": "Decreased", "correct": false},
        {"answer": "Stayed the same", "correct": false}
    ]}, {
    "question": "Compared to the partially simplified ARG, has the number of NODES in the fully simplified ARG increased, decreased, or stayed the same?",
    "type": "many_choice",
    "answers": [
        {"answer": "Increased", "correct": false},
        {"answer": "Decreased", "correct": false},
        {"answer": "Stayed the same", "correct": true, "feedback": "Correct: we have not deleted any nodes from the ARG, simply bypassed them if they appear unary in a local tree"}
    ]}, {
    "question": "Compared to the partially simplified ARG, has the number of MUTATIONS in the fully simplified ARG increased, decreased, or stayed the same??",
    "type": "many_choice",
    "answers": [
        {"answer": "Increased", "correct": false},
        {"answer": "Decreased", "correct": false},
        {"answer": "Stayed the same", "correct": true,  "feedback": "Correct, simplification doesn't change the pattern of encoded variation"}
    ]
}]

WB1_base["Q18.json"] = [{
    "question": "In general, when you go backwards in time, does the span of an ancestral haplotype in an ARG increase, decrease, or stay about the same?",
    "type": "many_choice",
    "answers": [
        {"answer": "Span decreases (gets shorter)", "correct": true, "feedback": "Correct: older haplotypes will have gone through more recombination events, and so the genomic region passed to the current-day samples will have been whittled down in length"},
        {"answer": "Span increases (gets longer)", "correct": false},
        {"answer": "Span stays about the same", "correct": false}
    ]
}]

WB1_base["Q19.json"] = [{
    "question": "How many sites now have no mutations associated with them?",
    "type": "numeric",
    "answers": [
        {"type": "value", "value": 1, "correct": true, "feedback": "Correct: the site at position 215 has no mutations."},
        {"type": "default"}
    ]
}]

WB1_base["Q20.json"] = [{
    "question": "Does this ARG with polytomies encode the same genetic data?",
    "type": "many_choice",
    "answers": [
        {"answer": "Yes: there is no change to the genetic data encoded", "correct": true},
        {"answer": "No", "correct": false},
        {"answer": "Sometimes", "correct": false},
    ]
}]

WB1_base["Q21.json"] = [{
    "question": "By looking at the provenance commands, which simulation software was used to make this ARG?",
    "type": "many_choice",
    "answers": [
        {"answer": "msprime", "correct": true},
        {"answer": "SLiM", "correct": false},
        {"answer": "fwdpy11", "correct": false},
        {"answer": "simuPOP", "correct": false},
        {"answer": "ms", "correct": false},
    ]
}]

WB1_base["Q22.json"] = [{
    "question": "From the log time plot, from when (in generations ago) does the exponential growth phase in both populations start model start, according to this model?",
    "type": "numeric",
    "answers": [
        {"type": "value", "value": 200, "correct": true},
        {"type": "default", "feedback": "Round to the nearest 100 generations"}
    ]
}]


WB1_base["Q25.json"] = [{
    "question": "Which of these random placement methods is most appropriate for choosing a node above which to place a neutral mutation?",
    "type": "many_choice",
    "answers": [
        {"answer": "With probability proportional to the span of the edge above the node", "correct": false},
        {"answer": "With equal probability above any node", "correct": false},
        {"answer": "Simply above the node with the longest branch", "correct": false},
        {"answer": "With probability proportional to the temporal length (evolutionary time) of the edge above the node", "correct": false},
        {"answer": "With probability proportional to the area (span × temporal length) of the edge above the node", "correct": true}
    ]
}]


def remove_edges(ts, edge_id_remove_list):
    edges_to_remove_by_child = collections.defaultdict(list)
    edge_id_remove_list = set(edge_id_remove_list)
    for m in ts.mutations():
        if m.edge in edge_id_remove_list:
            # If we remove this edge, we will remove the associated mutation
            # as the child node won't have ancestral material in this region.
            # So we force the user to explicitly (re)move the mutations beforehand
            raise ValueError("Cannot remove edges that have associated mutations")
    for remove_edge in edge_id_remove_list:
        e = ts.edge(remove_edge)
        edges_to_remove_by_child[e.child].append(e)

    # sort left-to-right for each child
    for k, v in edges_to_remove_by_child.items():
        edges_to_remove_by_child[k] = sorted(v, key=lambda e: e.left)
        # check no overlaps
        for e1, e2 in zip(edges_to_remove_by_child[k], edges_to_remove_by_child[k][1:]):
            assert e1.right <= e2.left

    # Sanity check: this means the topmost node will deal with modified edges left at the end
    assert ts.edge(-1).parent not in edges_to_remove_by_child
    
    new_edges = collections.defaultdict(list)
    tables = ts.dump_tables()
    tables.edges.clear()
    samples = set(ts.samples())
    # Edges are sorted by parent time, youngest first, so we can iterate over
    # nodes-as-parents visiting children before parents by using itertools.groupby
    for parent_id, ts_edges in itertools.groupby(ts.edges(), lambda e: e.parent):
        # Iterate through the ts edges *plus* the polytomy edges we created in previous steps.
        # This allows us to re-edit polytomy edges when the edges_to_remove are stacked
        edges = list(ts_edges)
        if parent_id in new_edges:
             edges += new_edges.pop(parent_id)
        if parent_id in edges_to_remove_by_child:
            for e in edges:
                assert parent_id == e.parent
                l = -1
                if e.id in edge_id_remove_list:
                    continue
                # NB: we go left to right along the target edges, reducing edge e as required
                for target_edge in edges_to_remove_by_child[parent_id]:
                    # As we go along the target_edges, gradually split e into chunks.
                    # If edge e is in the target_edge region, change the edge parent
                    assert target_edge.left > l
                    l = target_edge.left
                    if e.left >= target_edge.right:
                        # This target edge is entirely to the LHS of edge e, with no overlap
                        continue
                    elif e.right <= target_edge.left:
                        # This target edge is entirely to the RHS of edge e with no overlap.
                        # Since target edges are sorted by left coord, all other target edges
                        # are to RHS too, and we are finished dealing with edge e
                        tables.edges.append(e)
                        e = None
                        break
                    else:
                        # Edge e must overlap with current target edge somehow
                        if e.left < target_edge.left:
                            # Edge had region to LHS of target
                            # Add the left hand section (change the edge right coord)
                            tables.edges.add_row(left=e.left, right=target_edge.left, parent=e.parent, child=e.child)
                            e = e.replace(left=target_edge.left)
                        if e.right > target_edge.right:
                            # Edge continues after RHS of target
                            assert e.left < target_edge.right
                            new_edges[target_edge.parent].append(
                                e.replace(right=target_edge.right, parent=target_edge.parent)
                            )
                            e = e.replace(left=target_edge.right)
                        else:
                            # No more of edge e to RHS
                            assert e.left < e.right
                            new_edges[target_edge.parent].append(e.replace(parent=target_edge.parent))
                            e = None
                            break
                if e is not None:
                    # Need to add any remaining regions of edge back in 
                    tables.edges.append(e)
        else:
            # NB: sanity check at top means that the oldest node will have no edges above,
            # so the last iteration should hit this branch
            for e in edges:
                if e.id not in edge_id_remove_list:
                    tables.edges.append(e)
    assert len(new_edges) == 0
    tables.sort()
    return tables.tree_sequence()

def unsupported_edges(ts, per_interval=False):
    """
    Return the internal edges that are unsupported by a mutation.
    If ``per_interval`` is True, each interval needs to be supported,
    otherwise, a mutation on an edge (even if there are multiple intervals
    per edge) will result in all intervals on that edge being treated
    as supported.
    """
    edges_to_remove = np.ones(ts.num_edges, dtype="bool")
    edges_to_remove[[m.edge for m in ts.mutations()]] = False
    # We don't remove edges above samples
    edges_to_remove[np.isin(ts.edges_child, ts.samples())] = False

    if per_interval:
        return np.where(edges_to_remove)[0]
    else:
        keep = (edges_to_remove == False)
        for p, c in zip(ts.edges_parent[keep], ts.edges_child[keep]):
            edges_to_remove[np.logical_and(ts.edges_parent == p, ts.edges_child == c)] = False
        return np.where(edges_to_remove)[0]


class Workbook:
    @staticmethod
    def setup():
        display(HTML(
            "<style type='text/css'>" +
            ".exercise {background-color: yellow; color: black; font-family: 'serif'; font-size: 1.2em}" +
            ".exercise code {font-size: 0.7em}" +
            "</style>" + 
            "<h4>✅ Your notebook is ready to go!</h4>" +
            ("" if "pyodide" not in sys.modules else '''
To clear the notebook and reset,
select "Clear Browser Data" from the JupyterLite help menu.
''')
    ))

    def node_coalescence_status(arg):
        """
        Uses the num_children_array attribute to find nodes that represent local coalescence.
        See https://tskit.dev/tskit/docs/latest/python-api.html#tskit.Tree.num_children_array
        Returns an array of length num_nodes containing 0 if a node never has any coalescent
        segments, 1 if some segments of the node are coalescent and some unary, and 2 if
        all node segments represent a local coalescence point.
        """
        has_unary = np.zeros(arg.num_nodes + 1, dtype=int)
        has_coal = np.zeros(arg.num_nodes + 1, dtype=int)
        for tree in arg.trees():
            has_unary[tree.num_children_array == 1] = 1
            has_coal[tree.num_children_array > 1] = 1
        status = np.where(has_coal, np.where(has_unary, 1, 2), 0)
        return status[:-1] # remove the last array value, which is the "virtual root": see docs

    def remove_unsupported_edges(ts, per_interval=True):
        """
        Remove edges from the tree sequence that are unsupported by mutations.
        If ``per_interval`` is True, each interval needs to be supported,
        otherwise, a mutation on an edge (even if there are multiple intervals
        per edge) will result in all intervals on that edge being treated
        as supported.
        """
        edges_to_remove = unsupported_edges(ts, per_interval=per_interval)
        tables = remove_edges(ts, edges_to_remove).dump_tables()
        tables.edges.squash()
        return tables.tree_sequence()

class Workbook1(Workbook):
    url = WB1_base  # Put the real URL base string (ending in "/") here, once JSON question files have been made available at that URL
