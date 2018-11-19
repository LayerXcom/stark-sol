from src import fft, mimc_stark, merkle_tree, compression, fri
import json
# from mimc_stark import mk_mimc_proof, modulus, mimc, verify_mimc_proof
# from src.compression import compress_fri, compress_branches, bin_length
# from merkle_tree import merkelize, mk_branch, verify_branch
# from src.fri import prove_low_degree, verify_low_degree_proof

def test_merkletree():
    t = merkle_tree.merkelize([x.to_bytes(32, 'big') for x in range(128)])
    b = merkle_tree.mk_branch(t, 59)
    assert merkle_tree.verify_branch(t[1], 59, b, output_as_int=True) == 59
    print('Merkle tree works')
    
def test_fri():
    # Pure FRI tests
    poly = list(range(4096))
    root_of_unity = pow(7, (mimc_stark.modulus-1)//16384, mimc_stark.modulus)
    evaluations = fft.fft(poly, mimc_stark.modulus, root_of_unity)
    proof = fri.prove_low_degree(evaluations, root_of_unity, 4096, mimc_stark.modulus)
    print("Approx proof length: %d" % compression.bin_length(compression.compress_fri(proof)))
    assert fri.verify_low_degree_proof(merkle_tree.merkelize(evaluations)[1], root_of_unity, proof, 4096, mimc_stark.modulus)
    
    try:
        fakedata = [x if pow(3, i, 4096) > 400 else 39 for x, i in enumerate(evaluations)]
        proof2 = fri.prove_low_degree(fakedata, root_of_unity, 4096, mimc_stark.modulus)
        assert fri.verify_low_degree_proof(merkle_tree.merkelize(fakedata)[1], root_of_unity, proof, 4096, mimc_stark.modulus)
        raise Exception("Fake data passed FRI")
    except:
        pass
    try:
        assert fri.verify_low_degree_proof(merkle_tree.merkelize(evaluations)[1], root_of_unity, proof, 2048, mimc_stark.modulus)
        raise Exception("Fake data passed FRI")
    except:
        pass

def test_stark(testlang):
    INPUT = 3
    import sys
    # LOGSTEPS = int(sys.argv[1]) if len(sys.argv) > 1 else 13
    LOGSTEPS = 10
    # Full STARK test
    import random
    #constants = [random.randrange(mimc_stark.modulus) for i in range(64)]
    constants = [(i**7) ^ 42 for i in range(64)]
    proof = mimc_stark.mk_mimc_proof(INPUT, 2**LOGSTEPS, constants)
    m_root, l_root, branches, fri_proof = proof
    # L1 = compression.bin_length(compression.compress_branches(branches))
    # L2 = compression.bin_length(compression.compress_fri(fri_proof))
    # print("Approx proof length: %d (branches), %d (FRI proof), %d (total)" % (L1, L2, L1 + L2))
    # assert mimc_stark.verify_mimc_proof(3, 2**LOGSTEPS, constants, mimc_stark.mimc(3, 2**LOGSTEPS, constants), proof)    
    component0 = []
    component1 = []
    component2 = []
    for component in proof[3][:-1]:
        component0.append(component[0])    

    for component in proof[3][:-1]:
        for i in range(len(component[1])):
            component1.append(component[1][i])        

    for component in proof[3][:-1]:
        for i in range(len(component[2])):            
            for j in range(len(component[2][i])):
                component2.append(component[2][i][j])
    
    data = {
        'input': 3,
        'steps': 2**LOGSTEPS,
        'round_constants': constants,
        'output': mimc_stark.mimc(3, 2**LOGSTEPS, constants),
        'proof': {
            'root': proof[0].hex(),
            'lRoot': proof[1].hex(),
            'branches': [list(map(lambda x: x.hex(), branch)) for branch in proof[2]],
            'fri_components': {                                
                'root2': [list(map(lambda x: x.hex(), component0))],
                'branches2': {
                    'branchForColumns': [list(map(lambda x: x.hex(), componentForColumn)) for componentForColumn in component1],
                    'branchForPolys' : [list(map(lambda x: x.hex(), componentForPoly)) for componentForPoly in component2]
                },
                'directProof': [list(map(lambda x: x.hex(), proof[3][-1]))]
            }
        }
    }
    
    fw = open('stark_proof.json', 'w')
    json.dump(data, fw, indent='\t')
    assert testlang.verifier_contract.verifyMimcProof(3, 2**LOGSTEPS, constants, mimc_stark.mimc(3, 2**LOGSTEPS, constants), proof) == True

if __name__ == '__main__':
    test_stark()
