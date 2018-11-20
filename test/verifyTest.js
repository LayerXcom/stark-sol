
const proofData = require("./stark_proof.json");
const { input, steps, output, proof } = proofData;
const roundConstants = proofData.round_constants;

const BN = web3.utils.BN;

require('chai')
    .use(require('chai-as-promised'))    
    // .use(require('chai-bignumber')(BigNumber))
    .should();

let friComponents = [];

for (var i; i < proof.fri_components.root2.length; i++) {
    friComponents.push({
        root: proof.fri_components.root2[i],
        branchForColumns: proof.fri_components.branches2.branch_for_columns[i],
        branchesForPolys: proof.fri_components.branches2.branch_for_polys[i],
        directProof: ''
    })
}

friComponents.push({
    root: '',
    branchForColumns: '',
    branchesForPolys: '',
    directProof: proof.fri_components.direct_proof
})

let proofForStark = {
    root: proof.root,
    lRoot: proof.lRoot,
    branches: proof.branches,
    friComponent: friComponents
}

const Verifier = artifacts.require("VerifierContract");

contract('Verifier', () => {
    let verifier;

    beforeEach(async () => {
        verifier = await Verifier.new();
    });

    it("should be success", async () => {
        const isVerifiable = await verifier.verifyMimcProof(new BN(input), roundConstants, new BN(output), proofForStark); // TODO: 
        isVerifiable.should.equal(true);
    })
})