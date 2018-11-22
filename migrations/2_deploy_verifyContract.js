const Verifier = artifacts.require("VerifierContract");
const Merkle = artifacts.require("Merkle");

module.exports = async (deployer) => {    
    await deployer.deploy(Merkle)
    await deployer.link(Merkle, Verifier)
    await deployer.deploy(Verifier)                
        
}