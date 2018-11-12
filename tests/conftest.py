import os
import pytest
from ethereum import utils
from ethereum.tools import tester
from ethereum.abi import ContractTranslator
from ethereum.config import config_metropolis
from src.deployer import Deployer
from solc_simple import Builder

GAS_LIMIT = 8000000
START_GAS = GAS_LIMIT - 1000000
config_metropolis['BLOCK_GAS_LIMIT'] = GAS_LIMIT

# Compile contracts before testing
OWN_DIR = os.path.dirname(os.path.realpath(__file__))
CONTRACTS_DIR = os.path.abspath(os.path.realpath(os.path.join(OWN_DIR, '../contracts')))
OUTPUT_DIR = os.path.abspath(os.path.realpath(os.path.join(OWN_DIR, '../build')))
builder = Builder(CONTRACTS_DIR, OUTPUT_DIR)
builder.compile_all()
deployer = Deployer(builder)

@pytest.fixture
def ethtester():
    tester.chain = tester.Chain()
    return tester

@pytest.fixture
def ethutils():
    return utils

@pytest.fixture
def get_contract(ethtester, ethutils):
    def create_contract(path, args=(), sender=ethtester.k0):
        abi, hexcode = deployer.builder.get_contract_data(path)
        bytecode = ethutils.decode_hex(hexcode)
        encoded_args = (ContractTranslator(abi).encode_constructor_arguments(args) if args else b'')
        code = bytecode + encoded_args
        address = ethtester.chain.tx(sender=sender, to=b'', startgas=START_GAS, data=code)
        return ethtester.ABIContract(ethtester.chain, abi, address)
    return create_contract

@pytest.fixture
def verifier_contract(ethtester, get_contract):
    contract = get_contract('Verifier')
    ethtester.chain.mine()
    return contract