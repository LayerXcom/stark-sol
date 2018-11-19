
const proofData = require("./stark_proof.json");
const { input, steps, output, proof } = proofData;
const roundConstants = proofData.round_constants;


require('chai')
    .use(require('chai-as-promised'))    
    .should();

let friComponents = [];

for (var i; i<2; i++) {
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
    root: '',
    lRoot: '',
    branches: '',
    friComponent: [
        {
            root: '',
            branchForColumns: '',
            branchesForPolys: '',
            friComponent: ''
        }
    ]
}

proofForStark = {
    root: proof.root,
    lRoot: proof.lRoot,
    branches: proof.branches,
    friComponent: friComponents
}




const Verifier = artifacts.require("VerifierContract");

contract('Verifier', ([owner]) => {
    let verifier;

    beforeEach(async () => {
        verifier = await Verifier.new();
    });

    it("should be success", async () => {
        const isVerifiable = await verifier.verifyMimcProof(input, roundConstants, output, proofForStark); // TODO: 
        isVerifiable.should.equal(true);
    })

})