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
        bytes[][] branches;           // branches of P, D and B - evaluations
        FriComponent[] friComponent;  // low-degree proofs
    }

    struct FriComponent {
        bytes32 root;      // merkle root of columns
        // bytes[] branches;  // branches of the column and the four values in the polynominal
        bytes32[] branchForColumns;
        bytes32[] branchesForPolys;
        bytes32[] directProof;
    }
    
    struct Data {
        uint precision;
        uint G2;
        uint skips;
        uint skips2;
        uint lastStepPosition;
        uint[] constantsMiniPolynomial;
        bytes[][] branches;
        FriComponent[] friComponent;
        uint[] positions;
    }
    
    struct Data2 {
        uint x;
        uint xToThe_steps;
        bytes mBranch1;
        bytes mBranch2;
        uint lx;
        uint px;
        uint pG1x;
        uint dx;
        uint bx;
        uint zValue;
        uint kx;
        uint[] interpolant;
        uint[] zeropoly2;
    }

    // (for avoiding overflow) 2 ** 256 - 2 ** 32 * 351 + 1 =
    uint constant MODULUS = 115792089237316195423570985008687907853269984665640564039457584006405596119041;
    uint constant SPOT_CHECK_SECURITY_FACTOR = 80;
    uint constant EXTENSION_FACTOR = 8;

    uint _steps = 2 ** 3;
    uint G2 = (7 ** ((MODULUS - 1).div(_steps.mul(EXTENSION_FACTOR)))) % MODULUS;
    uint skips = _steps.mul(EXTENSION_FACTOR).div(_steps);
    uint lastStepPosition = (G2 ** ((_steps - 1).mul(skips))) % MODULUS;
    uint precision = _steps.mul(EXTENSION_FACTOR);

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
        // uint _steps, 
        uint[] _roundConstants, 
        uint _output, 
        Proof _proof        
    ) public returns (bool) 
    {
        require(_steps <= 2 ** 32);
        require(isPowerOf2(_steps) && isPowerOf2(_roundConstants.length));
        require(_roundConstants.length < _steps);
        
        Data memory data = Data({
            precision: _steps.mul(EXTENSION_FACTOR),
            G2:  (7 ** ((MODULUS - 1).div(_steps.mul(EXTENSION_FACTOR)))) % MODULUS,  // TODO: fix not to be overflowed
            skips: _steps.mul(EXTENSION_FACTOR).div(_steps),
            skips2: _steps.div(_roundConstants.length),
            lastStepPosition: (((7 ** ((MODULUS - 1).div(_steps.mul(EXTENSION_FACTOR)))) % MODULUS) ** ((_steps - 1).mul(_steps.div(_roundConstants.length)))) % MODULUS,
            constantsMiniPolynomial: fft(_roundConstants, MODULUS, (((7 ** ((MODULUS - 1).div(_steps.mul(EXTENSION_FACTOR)))) % MODULUS) ** (EXTENSION_FACTOR * (_steps.div(_roundConstants.length)))) % MODULUS, true),
            branches: _proof.branches,
            friComponent: _proof.friComponent,
            positions: getPseudorandomIndices(_proof.lRoot, _steps.mul(EXTENSION_FACTOR), SPOT_CHECK_SECURITY_FACTOR, EXTENSION_FACTOR)
        });
        
        require(verifyLowDegreeProof(_proof.lRoot, data.G2, data.friComponent, _steps * 2, MODULUS, EXTENSION_FACTOR));


        for (uint i; i < data.positions.length; i++) {
            
            Data2 memory data2 = Data2({
                x: (data.G2 ** data.positions[i]) % MODULUS,
                xToThe_steps: (((data.G2 ** data.positions[i]) % MODULUS) ** _steps) % MODULUS,
                mBranch1: _proof.root.verifyBranch(data.positions[i], data.branches[i * 3]),    // a branch check for P, D and B
                mBranch2: _proof.root.verifyBranch((data.positions[i].add(data.skips)) % _steps.mul(EXTENSION_FACTOR), data.branches[i * 3 + 1]),   // a branch check for P of g1x
                lx: _proof.root.verifyBranch(data.positions[i], data.branches[i * 3 + 2]).toUint(0),  // a branch check for L
                px: (_proof.root.verifyBranch(data.positions[i], data.branches[i * 3])).slice(0, 32).toUint(0),
                pG1x: (_proof.root.verifyBranch((data.positions[i].add(data.skips)) % _steps.mul(EXTENSION_FACTOR), data.branches[i * 3 + 1])).slice(0, 32).toUint(0),
                dx: (_proof.root.verifyBranch(data.positions[i], data.branches[i * 3])).slice(32, 32).toUint(0),
                bx: (_proof.root.verifyBranch((data.positions[i].add(data.skips)) % _steps.mul(EXTENSION_FACTOR), data.branches[i * 3 + 1])).slice(64, 32).toUint(0),
                zValue: polyDiv((((data.G2 ** data.positions[i]) % MODULUS) ** _steps) % MODULUS - 1, ((data.G2 ** data.positions[i]) % MODULUS) - data.lastStepPosition),
                kx: evalPolyAt(data.constantsMiniPolynomial, (((data.G2 ** data.positions[i]) % MODULUS) ** data.skips2) % MODULUS),
                interpolant: lagrangeInterp2([1, data.lastStepPosition], [_input, _output]),
                zeropoly2: mulPolys([uint(-1), 1], [-data.lastStepPosition, 1])
            });                        
    
            // Check transition constraints C(P(x)) = Z(x) * D(x)
            require((data2.pG1x - data2.px ** 3 - data2.kx - data2.zValue * data2.dx) % MODULUS == 0);
    
            // Check boundary constraints B(x) * Q(x) + I(x) = P(x)            
            require((data2.px - data2.bx * evalPolyAt(data2.zeropoly2, data2.x) - evalPolyAt(data2.interpolant, data2.x)) % MODULUS == 0);
    
            // Check correctness of the linear combination
            require((data2.lx - data2.dx - uint(keccak256(abi.encodePacked(_proof.root, 0x01))) * data2.px - uint(keccak256(abi.encodePacked(_proof.root, 0x02))) * data2.px * data2.xToThe_steps -  uint(keccak256(abi.encodePacked(_proof.root, 0x03))) * data2.bx - uint(keccak256(abi.encodePacked(_proof.root, 0x04))) * data2.bx * data2.xToThe_steps) % MODULUS == 0);
            return true;
        }

        return true;
    }
    
    function setNewStep(uint _newSteps) public {
        _steps = _newSteps;
    }

    function isPowerOf2(uint _x) public pure returns (bool) {
        if (_x % 2 == 0) {
            return isPowerOf2(_x.div(2));
        }

        if (_x == 1) {
            return true;
        }

        return false;
    }

    function getPseudorandomIndices(bytes32 _seed, uint _modulus, uint _count, uint _excludeMultiplesOf) internal returns (uint[]) {
        require(_modulus < 2**24);
        bytes memory data = abi.encodePacked(_seed);

        while (data.length < _count*4) {
            data.concat(abi.encodePacked(keccak256(data.slice(data.length - 32, 32))));
        }

        uint[] memory out = new uint[](data.length / 4);
        uint start = 0;
        uint j = 0;
        if (_excludeMultiplesOf == 0) {
            for (j = 0; j < out.length; j++) {
                out[j] = (data.slice(start, 4).toUint(0)).mod(_modulus);
                start += 4;
            }
        } else {
            uint realModulus = _modulus.mul(_excludeMultiplesOf - 1);
            for (j = 0; j < out.length; j++) {
                uint x = (data.slice(start, 4).toUint(0)).mod(realModulus);
                out[j] = x.add(1).add(x.div(_excludeMultiplesOf - 1));
                start += 4;
            }
        }
        return out;
    }

    function polyDiv(uint _x, uint _y) internal returns (uint) {
        return 0;
    }

    function evalPolyAt(uint[] _p, uint _x) internal returns (uint) {
        uint out = 0;
        uint powerOfX = 1;

        for (uint i = 0; i < _p.length; i++) {
            out = out.add(powerOfX.mul(_p[i]));
            powerOfX = powerOfX.mul(_x).mod(MODULUS);
        }
        return out.mod(MODULUS);
    }

    function _simple_ft(uint[] _vals, uint _modulus, uint[] _roots) internal returns (uint[]) {
        uint[] memory out = new uint[](_vals.length);
        uint L = _roots.length;

        for (uint i = 0; i < L; i++) {
            uint v = 0;
            for (uint j = 0; j < L; j++) {
                v = v.add(_vals[j].mul(_roots[(i.mul(j)).mod(L)]));
            }
            out[i] = v.mod(_modulus);
        }

        return out;
    }

    function _fft(uint[] _vals, uint _modulus, uint[] _roots) internal returns (uint[]) {
        if (_roots.length <= 4) {
            return _simple_ft(_vals, _modulus, _roots);
        }

        uint halfLength = _vals.length.div(2);

        uint[] memory L = new uint[](halfLength);
        uint[] memory R = new uint[](halfLength);
        uint[] memory _eRoots = new uint[](halfLength);

        for (uint i = 0; i < _vals.length; i++) {
            uint j = i.div(2);
            if (i.mod(2) == 0) {
                L[j] = _vals[i];
                _eRoots[j] = _roots[i];
            } else {
                R[j] = _vals[i];
            }
        }

        // recursive
        L = _fft(L, _modulus, _eRoots);
        R = _fft(R, _modulus, _eRoots);

        uint[] memory out = new uint[](_vals.length);

        for (i = 0; i < L.length; i++) {
            uint yTimesRoot = R[i].mul(_roots[i]);
            out[i] = (L[i].add(yTimesRoot)).mod(_modulus);
            out[i.add(L.length)] = (L[i].sub(yTimesRoot)).mod(_modulus);
        }
        return out;
    }

    function fft(uint[] _vals, uint _modulus, uint _rootOfUnity, bool isInv) internal returns (uint[]) {
        require(_vals.length.mod(2) == 0);

        // calculate the order
        uint order = 0;
        uint i = 1;
        while (i != 1) {
            order = order.add(1);
            i = i.mul(_rootOfUnity);
        }

        order = order.sub(1);

        // create an array of roots
        uint[] memory roots = new uint[](order);
        for (i = 0; i < roots.length; i++){
            uint root;
            if (i == 0) {
                root = 1;
            } else if (i == 1) {
                root = _rootOfUnity;
            } else {
                root = roots[i-1].mul(_rootOfUnity).mod(_modulus);
            }
            roots[i] = root;
        }

        uint[] memory arr = new uint[](roots.length);
        for (i = 0; i < roots.length; i++) {
            if (i < _vals.length) {
                arr[i] = _vals[i];
            } else {
                arr[i] = 0;
            }
        }

        _vals = arr;

        if (isInv) {
            for (i = 0; i < roots.length; i++) {
                arr[i] = roots[roots.length.sub(1).sub(i)];
            }

            uint[] memory xs = _fft(_vals, _modulus, arr);

            uint invLen = 1;
            for (i = 0; i < _modulus.sub(2); i++) {
                invLen = invLen.mul(invLen);
            }
            invLen = invLen.mod(_modulus);

            uint[] memory out = new uint[](xs.length);
            for (i = 0; i < xs.length; i++) {
                out[i] = (xs[i].mul(invLen)).mod(_modulus);
            }
            return out;
        }
        return _fft(_vals, _modulus, roots);
    }

    function lagrangeInterp(uint[] _xs, uint[] _ys) internal returns (uint[]) {
        uint[] memory a = new uint[](3);
        return a;
    }

    // optimized for degree 2
    function lagrangeInterp2(uint[2] _xs, uint[2] _ys) internal returns (uint[]) {
        uint[] memory a = new uint[](3);
        return a;
    }

    function multiInterp4(uint[][] _xsets, uint[][] _ysets) internal returns (uint[][]) {
        uint[][] memory a = new uint[][](3);
        return a;
    }

    function mulPolys(uint[2] _a, uint[2] _b) internal returns (uint[]) {
        uint[] memory out = new uint[]((_a.length).add(_b.length).sub(1));

        for (uint i = 0; i < _a.length; i++) {
            for (uint j = 0; j < _b.length; j++) {
                out[i+j] = _a[i] * _b[j];
            }
        }
        return out;
    }
}
