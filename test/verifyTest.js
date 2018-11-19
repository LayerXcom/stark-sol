
const proofData = require("./starl_proof.json");
const { input, steps, output, proof } = proofData;
const roundConstants = proofData.round_constants;

const BigNumber = web3.BigNumber;

require('chai')
    .use(require('chai-as-promised'))
    .use(require('chai-bignumber')(BigNumber))
    .should();

let fricomponents;

for (const i; i<2; i++) {
    fricomponents = [].push({
        root: proof.fri_component.root2[i],
        branchForColumns: proof.fri_component.branches2.branch_for_columns[i],
        branchesForPolys: proof.fri_component.branches2.branch_for_polys[i],
        directProof: ''
    })
}

fricomponents.push({
    root: proof.fri_component.direct_proof,
    branchForColumns: '',
    branchesForPolys: ''
})


const proof = {
    root: proof.root,
    lRoot: proof.lRoot,
    branches: proof.branches,
    friComponent: fricomponents
}


const Verifier = artifacts.require("Verifier");

contract('Verifier', ([owner]) => {

    let verifier;

    beforeEach(async () => {
        verifier = await Verifier.new();
    });

    it("should be success", async () => {
        const isVerifiable = await verifier.verifyMimcProof(input, roundConstants, output, proof);
        isVerifiable.should.equal(true);
    })

})