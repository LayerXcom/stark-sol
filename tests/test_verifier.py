import pytest

def test_is_power_of_2(testlang):
    assert testlang.verifier_contract.isPowerOf2(16) == True    