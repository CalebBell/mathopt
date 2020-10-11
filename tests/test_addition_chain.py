from fluids.numerics import assert_close, assert_close1d
from mathopt.addition_chain import *

def test_tabulated_addition_chains():
    assert 1 not in tabulated_addition_chains
    assert 0 not in tabulated_addition_chains
    assert [[2, 4]] == tabulated_addition_chains[4]
    assert [[2, 3, 6, 12], [2, 4, 6, 12], [2, 4, 8, 12]] == tabulated_addition_chains[12]
    assert len(tabulated_addition_chains) > 500
    
def test_addition_chain_length():
    assert 8 == addition_chain_length(128)


def test_minimum_addition_chain_multi():
    calc = minimum_addition_chain_multi([3, 5, 11])
#    expect = (5, [[2, 3], [2, 3, 5], [2, 3, 5, 10, 11]])
#    assert calc == expect
    assert calc[0] == 5
    
    calc = minimum_addition_chain_multi([2, 3, 9, 29])
#    assert calc == (7, [[2], [2, 3], [2, 3, 6, 9], [2, 3, 6, 9, 18, 27, 29]])
    assert calc[0] == 7

def test_minimum_addition_chain_multi_heuristic():
    powers = [3, 5, 11, 12, 30, 100]
    l0, calc = minimum_addition_chain_multi_heuristic(powers, small_chain_length=5)
    l1, _ =  minimum_addition_chain_multi(powers)
    assert l0 == l1
    assert calc == [[2, 3],
                      [2, 3, 5],
                      [2, 3, 5, 6, 11],
                      [2, 3, 6, 12],
                      [2, 3, 6, 12, 24, 30],
                      [2, 3, 6, 12, 24, 48, 50, 100]]