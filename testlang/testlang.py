from ethereum.utils import sha3
import conftest


def get_accounts(ethtester):
    """Converts ethereum.tools.tester accounts into a list.
    Args:
        ethtester (ethereum.tools.tester): Ethereum tester instance.
    Returns:
        EthereumAccount[]: A list of EthereumAccounts.
    """

    accounts = []
    for i in range(10):
        address = getattr(ethtester, 'a{0}'.format(i))
        key = getattr(ethtester, 'k{0}'.format(i))
        accounts.append(EthereumAccount(address_to_hex(address), key))
    return accounts


class TestingLanguage(object):

    def __init__(self, verifier_contract, ethtester):
        self.verifier_contract = verifier_contract
        self.ethtester = ethtester
        self.accounts = get_accounts(ethtester)
