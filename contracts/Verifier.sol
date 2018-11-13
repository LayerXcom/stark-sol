pragma solidity ^0.4.24;
pragma experimental ABIEncoderV2;

contract VerifierContract {

    struct Proof {
        bytes32 root;               // merkle root of P, D and B - evaluations 
        bytes32 lRoot;              // merkle root of L - evaluations
        bytes[] branches;           // branches of P, D and B - evaluations
        FriComponent friComponent;  // low-degree proofs
    }

    struct FriComponent {
        bytes32 root;      // merkle root of columns
        bytes[] branches;  // branches of the column and the four values in the polynominal
    }


    // verify an FRI proof
    function verifyLowDegreeProof(
        bytes32 merkleRoot, 
        uint rootOfUnity, 
        FriComponent[] friComponents, 
        uint maxDegPlus1, 
        uint modulus, 
        uint excludeMultiplesOf
    ) public returns (bool) {
        return false;
    }

    // verify a STARK
    function verifyMimcProof(
        uint input, 
        uint steps, 
        uint roundConstants, 
        uint output, 
        Proof proof
    ) public returns (bool) {
        return false;
    }
}
