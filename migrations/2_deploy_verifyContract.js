const Verifier = artifacts.require("VerifierContract");

module.exports = (deployer) => {
    deployer.deploy(Verifier);
}