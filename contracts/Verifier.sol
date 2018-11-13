pragma solidity ^0.4.24;
pragma experimental ABIEncoderV2;

import "./SafeMath.sol";
import "./BytesLib.sol";
import "./Merkle.sol";

contract VerifierContract {
    using SafeMath for uint;
    using BytesLib for bytes;
    using Merkle for bytes32;

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

    uint constant MODULUS = 2 ** 256 - 2 ** 32 * 351 + 1;
    uint constant SPOT_CHECK_SECURITY_FACTOR = 80;
    uint constant EXTENSION_FACTOR = 8;


    // verify an FRI proof
    function verifyLowDegreeProof(
        bytes32 _merkleRoot, 
        uint _rootOfUnity, 
        FriComponent[] _friComponents, 
        uint _maxDegPlus1, 
        uint _modulus, 
        uint _excludeMultiplesOf
    ) public returns (bool) 
    {
        return true;
    }

    // verify a STARK
    function verifyMimcProof(
        uint _input, 
        uint _steps, 
        uint[] _roundConstants, 
        uint _output, 
        Proof _proof
    ) public returns (bool) 
    {
        bytes32 root = _proof.root;
        bytes32 lRoot = _proof.lRoot;
        bytes[] memory branches = _proof.branches;
        FriComponent memory friComponent = _proof.friComponent;

        require(_steps <= 2 ** 32);
        require(isPowerOf2(steps) && isPowerOf2(_roundConstants.length));
        require(_roundConstants.length < _steps);

        uint precision = _steps.mul(EXTENSION_FACTOR);
        uint G2 = (7 ** ((MODULUS - 1).div(precision))) % MODULUS;
        uint skips = precision.div(_steps);
        uint skips2 = _steps.div(_roundConstants.length);
        
        uint[] constantsMiniPolynomial = fft(_roundConstants, MODULUS, (G2 ** (EXTENSION_FACTOR * skips2) % MODULUS, true);
        require(verifyLowDegreeProof(lRoot, G2, friComponent, _steps * 2, MODULUS, EXTENSION_FACTOR));

        // uint k1 = keccak256(abi.encodePacked(root, 0x01)).toUint(0);
        // uint k2 = keccak256(abi.encodePacked(root, 0x02)).toUint(0);
        // uint k3 = keccak256(abi.encodePacked(root, 0x03)).toUint(0);
        // uint k4 = keccak256(abi.encodePacked(root, 0x04)).toUint(0);

        uint k1 = uint(keccak256(abi.encodePacked(root, 0x01)));
        uint k2 = uint(keccak256(abi.encodePacked(root, 0x02)));
        uint k3 = uint(keccak256(abi.encodePacked(root, 0x03)));
        uint k4 = uint(keccak256(abi.encodePacked(root, 0x04)));

        uint[] positions = getPseudoramdomIndicies(lRoot, precision, SPOT_CHECK_SECURITY_FACTOR, EXTENSION_FACTOR);
        uint lastStepPosition = (G2 ** ((steps - 1).mul(skips))) % MODULUS;

        for (uint i; i < positions.length; i++) {
            uint x = (G2 ** positions[i]) % MODULUS;
            uint xToTheSteps = (x ** _steps) % MODULUS;

            // a branch check for P, D and B
            bytes memory mBranch1 = root.verifyBranch(positions[i], branches[i * 3]);
            // a branch check for P of g1x
            bytes memory mBranch2 = root.verifyBranch((positions[i].add(skips)) % precision, branches[i * 3 + 1]);
            // a branch check for L
            uint lx = uint(root.verifyBranch(positions[i], branches[i * 3 + 2]));

            uint px = mBranch1.slice(0, 32).toUint(0);
            uint pG1x = mBranch2.slice(0, 32).toUint(0);
            uint dx = mBranch1.slice(32, 32).toUint(0);
            uint bx = mBranch2.slice(64, 32).toUint(0);

            uint zValue = polyDiv((x ** steps) % MODULUS - 1, x - lastStepPosition);
            uint kx = evalPolyAt(constantsMiniPolynomial, (x ** skips2) % MODULUS);

            // Check transition constraints C(P(x)) = Z(x) * D(x)
            require((pG1x - px ** 3 - kx - zValue * dx) % MODULUS == 0);

            // Check boundary constraints B(x) * Q(x) + I(x) = P(x)
            uint[3] interpolant = lagrangeInterp2([1, lastStepPosition], [_input, _output]);
            uint[] zeropoly2 = mulPolys([-1, 1], [-lastStepPosition, 1]);
            require((px - bx * evalPolyAt(zeropoly2, x) - evalPolyAt(interpolant, x)) % MODULUS == 0);

            // Check correctness of the linear combination
            require((lx - dx - k1 * px - k2 * px * xToTheSteps - k3 * bx - k4 * bx * xToTheSteps) % MODULUS == 0);
        }

        return true;
    }

    function isPowerOf2(uint _x) internal pure returns (uint) {
        return 0;
    }

    function getPseudorandomIndices(bytes32 _seed, uint _modulus, uint _count, uint _excludeMultiplesOf) internal returns (uint[]) {
        return [0];
    }

    function polyDiv(uint _x, uint _y) internal returns (uint) {
        return 0;
    }

    function evalPolyAt(uint[] _p, uint _x) internal returns (uint) {
        return 0;
    }

    function fft(uint[] _vals, uint _modulus, uint _rootOfUnity, bool isInv) internal returns (uint[]) {
        return [0];
    }

    function lagrangeInterp2(uint[2] _xs, uint[2] _ys) internal returns (uint[]) {
        return [0];
    }

    function mulPolys(uint[] _a, uint[] _b) internal returns (uint[]) {
        return [0];
    }

}
