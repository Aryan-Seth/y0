import sys

from collections.abc import Sequence
from utils import ZIdentification,Identification, Unidentifiable,Z2Identification
from y0.dsl import Expression, P, Probability, Product, Sum, Variable
from y0.graph import NxMixedGraph
from copy import deepcopy
__all__ = [
    "gz-identify",
]


def line_1(identification: ZIdentification):
    """
    Sum of probabilities of all vertices except target
    """
    outcomes = identification.outcomes
    vertices = set(identification.graph.nodes())
    return Sum.safe(
        expression=identification.estimand,
        ranges=vertices.difference(outcomes),
    )
    
def line_4(identification: ZIdentification) -> list[ZIdentification]:
    r"""Run line 4 of the identification algorithm.

    The key line of the algorithm, it decomposes the problem into a set
    of smaller problems using the key property of *c-component
    factorization* of causal models. If the entire graph is a single
    C-component already, further problem decomposition is impossible,
    and we must provide base cases. :math:`\mathbf{ID}` has three base
    cases.

    :param identification: The data structure with the treatment, outcomes, estimand, and graph
    :returns: A list of new estimands
    :raises ValueError: If the precondition that there are more than 1 districts without treatments is not met
    """
    treatments = identification.treatments
    estimand = identification.estimand
    graph = identification.graph
    vertices = set(graph.nodes())
    Z=identification.Z
    # line 4
    graph_without_treatments = graph.remove_nodes_from(treatments.union(identification.queryI.treatments).union(identification.queryJ.treatments))
    districts_without_treatment = graph_without_treatments.districts()
    if len(districts_without_treatment) <= 1:
        raise ValueError("Line 4 precondition not met")
    return [
        ZIdentification.from_parts(
            outcomes=set(district_without_treatment),
            treatments=vertices - district_without_treatment-Z,
            treatmentsI=identification.queryI.treatments,
            treatmentsJ=identification.queryJ.treatments.union(Z.intersection(vertices-district_without_treatment)),
            Z=Z-vertices-district_without_treatment,
            estimand=estimand,
            graph=graph,
        )
        for district_without_treatment in districts_without_treatment
    ]


def _get_single_district(graph: NxMixedGraph) -> frozenset[Variable]:
    districts = graph.districts()
    if len(districts) != 1:
        raise RuntimeError
    return districts.pop()

def line_7(identification: Identification) -> Identification:
    r"""Run line 7 of the identification algorithm.
    """
    outcomes = identification.outcomes
    treatments = identification.treatments
    graph = identification.graph

    graph_without_treatments = graph.remove_nodes_from(treatments)
    # line 7 precondition requires single district
    district_without_treatments = _get_single_district(graph_without_treatments)

    # line 7
    for district in graph.districts():
        if district_without_treatments < district:
            parents = list(graph.topological_sort())
            return Identification.from_parts(
                outcomes=outcomes,
                treatments=treatments.intersection(district),
                estimand=Product.safe(p_parents(v, parents) for v in district),
                graph=graph.subgraph(district),
            )

    raise ValueError("Could not identify suitable district")



def p_parents(child: Variable, ordering: Sequence[Variable]) -> Probability:
    """Get a probability expression based on a topological ordering.

    :param child: The child variable
    :param ordering: A topologically ordered sequence of all variables. All occurring before the
        child will be used as parents.
    :return: A probability expression
    """
    return P(child | ordering[: ordering.index(child)])

def gz_identify(identification: ZIdentification):
    """
    Run the gz-identification algorithm from z-transportability, Lee and Honavar

    """

    graph=identification.graph
    treatments = identification.treatments
    I= identification.treatmentsI
    J= identification.treatmentsJ
    Z= identification.Z
    outcomes = identification.outcomes
    vertices = set(graph.nodes())    

    #line 1
    if (len(treatments)==0 or treatments is None):
        return line_1(identification) #sum of probabilities of all vertices except target
    
    #line 2
    if(len(vertices-graph.ancestors_inclusive(outcomes))!=0):

        return gz_identify(ZIdentification.from_parts(
            outcomes=outcomes,
            treatments=treatments.intersection(graph.ancestors_inclusive(outcomes)),
            Z=Z,
            treatmentsI=I,
            treatmentsJ=J,
            graph=graph.subgraph(graph.ancestors_inclusive(outcomes)),
        ))
    
    #line 3
    W= set(vertices-(treatments.union(I).union(J))-graph.get_intervened_ancestors(treatments.union(I).union(J),outcomes=outcomes))
    Z_w= Z.intersection(W.union(treatments))

    if(Z_w.union(W) is None or len(Z_w.union(W))==0):
        return gz_identify(ZIdentification.from_parts(
            graph=graph,
            outcomes=outcomes,
            treatments=treatments.union(W)-Z_w,
            Z=Z - W,
            treatmentsI=I.union(Z_w),
            treatmentsJ=J,
        ))
    
    #line 4
    graph_without_treatments = graph.remove_nodes_from(treatments.union(I).union(J))
    ccomp_temp=graph_without_treatments.districts()

    if(len(ccomp_temp)>1):
        to_multiply=Product.safe(map(gz_identify,line_4(identification)))
        return Sum.safe(expression=to_multiply,ranges=vertices-treatments.union(I).union(outcomes))


    #line 5
    graph_without_treatments = graph.remove_nodes_from(treatments.union(I).union(J))
    if(len(set(graph_without_treatments.districts()))==1):
        return Unidentifiable(graph.nodes())
    
    #line 6
    if(graph_without_treatments.districts() in graph.districts() and len(graph_without_treatments.districts())==1):
        district=graph_without_treatments.districts().pop()
        parents=list(graph.topological_sort())
        expression=Product.safe(
            p_parents(v,parents-I.union(J)) for v in district
        )
        ranges=district-outcomes
        return Sum.safe(expression,ranges)
    
    #line 7
    return gz_identify(line_7(identification))
    

    

    
    