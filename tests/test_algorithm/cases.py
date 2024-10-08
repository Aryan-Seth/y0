"""Test cases."""

import unittest
from collections import Counter
from collections.abc import Collection

from y0.dsl import Expression, Variable, get_outcomes_and_treatments
from y0.graph import NxMixedGraph
from y0.mutate import canonicalize

__all__ = ["GraphTestCase"]


class GraphTestCase(unittest.TestCase):
    """Tests parallel worlds and counterfactual graphs."""

    def assert_graph_equal(
        self, expected: NxMixedGraph, actual: NxMixedGraph, msg=None, *, sort: bool = False
    ) -> None:
        """Check the graphs are equal (more nice than the builtin :meth:`NxMixedGraph.__eq__` for testing)."""
        if sort:
            self.assertEqual(
                sorted(set(expected.directed.nodes())),
                sorted(set(actual.directed.nodes())),
                msg=msg,
            )
        else:
            self.assertEqual(set(expected.directed.nodes()), set(actual.directed.nodes()), msg=msg)
        self.assertEqual(set(expected.undirected.nodes()), set(actual.undirected.nodes()), msg=msg)
        self.assertEqual(set(expected.directed.edges()), set(actual.directed.edges()), msg=msg)
        self.assertEqual(
            set(map(frozenset, expected.undirected.edges())),
            set(map(frozenset, actual.undirected.edges())),
            msg=msg,
        )

    def assert_expr_equal(self, expected: Expression, actual: Expression) -> None:
        """Assert that two expressions are the same."""
        expected_outcomes, expected_treatments = get_outcomes_and_treatments(query=expected)
        actual_outcomes, actual_treatments = get_outcomes_and_treatments(query=actual)
        self.assertEqual(expected_treatments, actual_treatments)
        self.assertEqual(expected_outcomes, actual_outcomes)
        ordering = sorted(expected.get_variables(), key=lambda x: str(x))
        expected_canonical = canonicalize(expected, ordering)
        actual_canonical = canonicalize(actual, ordering)
        self.assertEqual(
            expected_canonical,
            actual_canonical,
            msg=f"\nExpected: {expected_canonical!s}\nActual:   {actual_canonical!s}",
        )

    def assert_collection_of_set_equal(
        self, left: Collection[set[Variable]], right: Collection[set[Variable]]
    ) -> None:
        """Check that two collections contain sets with the same elements."""
        c1 = Counter(frozenset(element) for element in left)
        c2 = Counter(frozenset(el) for el in right)
        self.assertEqual(c1, c2)
