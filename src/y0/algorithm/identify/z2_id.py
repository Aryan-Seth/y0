"""An implementation of the identification algorithm."""
import sys

from collections.abc import Sequence
from utils import ZIdentification,Identification, Unidentifiable, Z2Identification
from y0.dsl import Expression, P, Probability, Product, Sum, Variable
from y0.graph import NxMixedGraph
from id_std import identify
from copy import deepcopy
__all__ = [
    "z2id",'subz2'
]

def subz2(identification: Z2Identification)->Z2Identification:
    graph=identification.graph
    treatment=identification.treatments
    outcome=identification.outcomes
    set_treatment = identification.set_treatments
    vertices = set(identification.graph.nodes())
    #line 1
    for i in range(len(set_treatment)):
        z= set(set_treatment[i])
        c_components=graph.districts()

        #condition2
        if(z<=treatment):
            for j in range(len(set_treatment)):
                z1=set(set_treatment[j])
                if(z<z1 and z1<=treatment):
                    return Unidentifiable(graph.nodes())
        return identify(identification)



    #line 2
    if(len(vertices-graph.ancestors_inclusive(outcome))!=0):
        outcomes_and_ancestors = graph.ancestors_inclusive(outcome)
        out_anc_graph=graph.subgraph(outcomes_and_ancestors)
        return subz2(Z2Identification.from_parts(
            outcomes=outcome,
            treatments=treatment.intersection(outcomes_and_ancestors),
            graph=out_anc_graph,
            set_treatments=set_treatment
        ))
    
    #line 3
    W=vertices-treatment-graph.get_intervened_ancestors(treatment,outcomes=outcome)
    if(len(W)!=0):
       return subz2(Z2Identification.from_parts(
         graph=graph,
         outcomes=outcome,
         treatments=treatment.union(W),
         set_treatments=set_treatment  
       ))
    
    #line 4
    subgraph=graph.remove_nodes_from(treatment)
    c_components=subgraph.districts()
    if(len(c_components)>1):
        to_multiply=[Z2Identification.from_parts(
            outcomes=component,
            treatments=vertices-set(component),
            graph=graph,
            set_treatments=set_treatment
        ) for component in c_components]
        expression=Product.safe(map(subz2,to_multiply))
        return Sum.safe(expression=expression,
                    ranges=vertices-treatment.union(outcome))
    return Unidentifiable(graph.nodes())

        

    
       
        
        
    

def z2id(identification: Z2Identification)->Z2Identification:
    graph=identification.graph
    treatment=identification.treatments
    outcome=identification.outcomes
    set_treatment=identification.set_treatments

    
    if(len(set_treatment[0])!=0 and len(treatment)==0):
        c_components=graph.districts()
        to_prod=[Z2Identification.from_parts(
            outcomes=district,
            treatments=set(graph.nodes())-set(district),
            graph=graph,
            set_treatments=set_treatment
        ) for district in c_components]

        expression=Product.safe(map(subz2,to_prod))
        return Sum.safe(expression=expression,ranges=set(identification.graph.nodes())-identification.outcomes)
    else:
        return subz2(identification)



