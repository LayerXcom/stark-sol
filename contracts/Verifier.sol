// pragma solidity ^0.4.24;
pragma solidity >=0.4.24 <0.6.0;
pragma experimental ABIEncoderV2;

import "./SafeMath.sol";
import "./BytesLib.sol";
import "./Merkle.sol";

contract VerifierContract {
    using SafeMath for uint;
    using BytesLib for bytes;
    using Merkle for bytes32;
    using Merkle for uint[];

    struct Proof {
        bytes32 root;               // merkle root of P, D and B - evaluations 
        bytes32 lRoot;              // merkle root of L - evaluations
        bytes[][] branches;           // branches of P, D and B - evaluations
        FriComponent[] friComponent;  // low-degree proofs
    }

    struct FriComponent {
        bytes32 root;      // merkle root of columns
        bytes[] branchForColumns;
        bytes[] branchesForPolys;
        bytes[] directProof;
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
        uint xToSteps;
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

    struct Data3 {
        bytes32 merkleRoot;
        uint rootOfUnity;
        uint maxDegPlus1;
        uint roudeg;
        uint[4] quarticRootsOfUnity;        
        uint excludeMultiplesOf;        
    }

    struct Data4 {
        uint[] ys;
        uint[4][] xcoords;
        uint[4][] rows;
        uint[] columnvals;
        uint j;
    }    

    // (for avoiding overflow) 2 ** 256 - 2 ** 32 * 351 + 1 =
    uint constant MODULUS = 115792089237316195423570985008687907853269984665640564039457584006405596119041;
    uint constant SPOT_CHECK_SECURITY_FACTOR = 80;
    uint constant EXTENSION_FACTOR = 8;

    // (7 ** (MODULUS / 16384)) % MODULUS    
    uint constant ROOTOFUNITY = 109030522278086733194866766814799573444537703137521227146756705386680725962607;

    uint _steps = 2 ** 3;
    uint G2 = (7 ** ((MODULUS - 1).div(_steps.mul(EXTENSION_FACTOR)))) % MODULUS;
    uint skips = _steps.mul(EXTENSION_FACTOR).div(_steps);
    uint lastStepPosition = (G2 ** ((_steps - 1).mul(skips))) % MODULUS;
    uint precision = _steps.mul(EXTENSION_FACTOR);    
    

        // verify an FRI proof
    function verifyLowDegreeProof(
        bytes32 _merkleRoot,  
        uint _rootOfUnity,
        FriComponent[] memory _friComponents,      
        uint _maxDegPlus1,         
        uint _excludeMultiplesOf
    ) internal returns (bool) 
    {          
        Data3 memory data3 = Data3({
            merkleRoot: _merkleRoot,
            rootOfUnity: _rootOfUnity,
            maxDegPlus1: _maxDegPlus1,
            roudeg: getRoudeg(),
            // Powers of the given root of unity 1, p, p**2, p**3 such that p**4 = 1
            quarticRootsOfUnity: [1, (_rootOfUnity ** (getRoudeg().div(4))) % MODULUS, (_rootOfUnity ** (getRoudeg().div(2))) % MODULUS, (_rootOfUnity ** (getRoudeg().mul(3).div(4))) % MODULUS],                                                                                                                        
            excludeMultiplesOf: _excludeMultiplesOf
        });

        for (uint i = 0; i < _friComponents.length - 1; i++) {
            data3.merkleRoot = verifyFriProofs(i, _friComponents, data3);
            data3.rootOfUnity = (data3.rootOfUnity ** 4) % MODULUS;
            data3.maxDegPlus1 = data3.maxDegPlus1.div(4);
            data3.roudeg = data3.roudeg.div(4);
        }                

        require(verifyDirectFriProof(_maxDegPlus1, _friComponents[_friComponents.length - 1].directProof, _merkleRoot, _excludeMultiplesOf));

        return true;
    }        

    function verifyFriProofs(uint i, FriComponent[] memory _friComponents, Data3 memory _data3) private returns (bytes32) {
        Data4 memory data4 = Data4({
            ys: getPseudorandomIndices(_friComponents[i].root, _data3.roudeg.div(4), 10, _data3.excludeMultiplesOf),
            xcoords: new uint[4][](getPseudorandomIndices(_friComponents[i].root, _data3.roudeg.div(4), 10, _data3.excludeMultiplesOf).length),
            rows: new uint[4][](getPseudorandomIndices(_friComponents[i].root, _data3.roudeg.div(4), 10, _data3.excludeMultiplesOf).length),
            columnvals: new uint[](getPseudorandomIndices(_friComponents[i].root, _data3.roudeg.div(4), 10, _data3.excludeMultiplesOf).length),
            j: 0
        });

        for (data4.j = 0; data4.j < data4.ys.length; data4.j++) {
            uint x1 = (_data3.rootOfUnity ** data4.ys[data4.j]) % MODULUS;
            uint[4] memory xcoord;
            uint[4] memory row;

            for (uint k = 0; k < 4; k++) {                                                    
                xcoord[k] = _data3.quarticRootsOfUnity[k].mul(x1) % MODULUS;                    
                row[k] = _data3.merkleRoot.verifyBranch(data4.ys[data4.j] + _data3.roudeg.div(4) * k, _friComponents[i].branchesForPolys).toUint(0); // TODO: devide 4 elements 
            }                

            data4.columnvals[data4.j] = _friComponents[i].root.verifyBranch(data4.ys[data4.j], _friComponents[i].branchForColumns).toUint(0); // TODO: devide 4 elements
            data4.xcoords[data4.j] = xcoord;
            data4.rows[data4.j] = row;
        }
        
        uint[][] memory polys = multiInterp4(data4.xcoords, data4.rows);
            
        for (data4.j = 0; data4.j < polys.length; data4.j++) { 
            // evaluate each polynomials at special x coord
            require(evalPolyAt(polys[data4.j], abi.encodePacked(_data3.merkleRoot).toUint(0) % MODULUS) == data4.columnvals[data4.j]);
        }
        return _friComponents[i].root;
    }

    function verifyDirectFriProof(uint _maxDegPlus1, bytes[] memory _directProof, bytes32 _merkleRoot, uint _excludeMultiplesOf) private returns (bool) {
        uint[] memory data = new uint[](_directProof.length);
        uint i;        
        for (i = 0; i < _directProof.length; i++) {
            data[i] = (_directProof[i]).toUint(0);            
        }

        require(_maxDegPlus1 <= 16);
        
        bytes32[] memory mtree = data.merkelize();
        require(mtree[1] == _merkleRoot);

        uint[] memory powers = getPowerCycle(ROOTOFUNITY);
        uint[] memory pts = new uint[](data.length);

        if (_excludeMultiplesOf != 0) {            
            for (i = 0; i < data.length; i++) {
                if (i % _excludeMultiplesOf != 0) {
                    pts[i] = i;
                }
            }
        } else {
            for (i = 0; i < data.length; i++) {
                pts[i] = i;
            }
        }

        uint[] memory slicedPts = new uint[](_maxDegPlus1);
        
        for (i = 0; i < _maxDegPlus1; i++) {
            slicedPts[i] = pts[i];            
        }

        uint[] memory xs2 = new uint[](slicedPts.length);
        uint[] memory ys2 = new uint[](slicedPts.length);
        
        for (i = 0; i < slicedPts.length; i++) {
            xs2[i] = powers[i];            
            ys2[i] = data[i];            
        }
        
        uint[] memory poly = lagrangeInterp(xs2, ys2);
        
        for (i = 0; i < slicedPts.length; i++) {
            require(evalPolyAt(poly, powers[i]) == data[i]);
        }

        return true;
    }

    // verify a STARK
    function verifyMimcProof(
        uint _input, 
        // uint _steps, 
        uint[] memory _roundConstants, 
        uint _output, 
        Proof memory _proof        
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
        
        require(verifyLowDegreeProof(_proof.lRoot, G2, data.friComponent, _steps * 2, EXTENSION_FACTOR));


        for (uint i; i < data.positions.length; i++) {
            
            Data2 memory data2 = Data2({
                x: (data.G2 ** data.positions[i]) % MODULUS,
                xToSteps: (((data.G2 ** data.positions[i]) % MODULUS) ** _steps) % MODULUS,
                mBranch1: _proof.root.verifyBranch(data.positions[i], data.branches[i * 3]),    // a branch check for P, D and B
                mBranch2: _proof.root.verifyBranch((data.positions[i].add(data.skips)) % _steps.mul(EXTENSION_FACTOR), data.branches[i * 3 + 1]),   // a branch check for P of g1x
                lx: _proof.root.verifyBranch(data.positions[i], data.branches[i * 3 + 2]).toUint(0),  // a branch check for L
                px: (_proof.root.verifyBranch(data.positions[i], data.branches[i * 3])).slice(0, 32).toUint(0),
                pG1x: (_proof.root.verifyBranch((data.positions[i].add(data.skips)) % _steps.mul(EXTENSION_FACTOR), data.branches[i * 3 + 1])).slice(0, 32).toUint(0),
                dx: (_proof.root.verifyBranch(data.positions[i], data.branches[i * 3])).slice(32, 32).toUint(0),
                bx: (_proof.root.verifyBranch((data.positions[i].add(data.skips)) % _steps.mul(EXTENSION_FACTOR), data.branches[i * 3 + 1])).slice(64, 32).toUint(0),
                zValue: div((((data.G2 ** data.positions[i]) % MODULUS) ** _steps) % MODULUS - 1, ((data.G2 ** data.positions[i]) % MODULUS) - data.lastStepPosition),
                kx: evalPolyAt(data.constantsMiniPolynomial, (((data.G2 ** data.positions[i]) % MODULUS) ** data.skips2) % MODULUS),
                interpolant: lagrangeInterp2([1, data.lastStepPosition], [_input, _output]),
                zeropoly2: mulPolys([uint(-1), 1], [-data.lastStepPosition, 1])
            });                        
    
            // Check transition constraints C(P(x)) = Z(x) * D(x)
            require((data2.pG1x - data2.px ** 3 - data2.kx - data2.zValue * data2.dx) % MODULUS == 0);
    
            // Check boundary constraints B(x) * Q(x) + I(x) = P(x)            
            require((data2.px - data2.bx * evalPolyAt(data2.zeropoly2, data2.x) - evalPolyAt(data2.interpolant, data2.x)) % MODULUS == 0);
    
            // Check correctness of the linear combination
            require((data2.lx - data2.dx - uint(keccak256(abi.encodePacked(_proof.root, "0x01"))) * data2.px - uint(keccak256(abi.encodePacked(_proof.root, "0x02"))) * data2.px * data2.xToSteps -  uint(keccak256(abi.encodePacked(_proof.root, "0x03"))) * data2.bx - uint(keccak256(abi.encodePacked(_proof.root, "0x04"))) * data2.bx * data2.xToSteps) % MODULUS == 0);
            return true;
        }

        return true;
    }
    
    function setNewStep(uint _newSteps) internal {
        _steps = _newSteps;
    }

    function getRoudeg() internal pure returns (uint) {
        uint testVal = ROOTOFUNITY;
        uint roudeg = 1;
        while( testVal != 1) {
            roudeg *= 2;
            testVal = (testVal * testVal) % MODULUS;
        }
        return roudeg;
    }

    uint[] o; // TODO: replace to memory val

    function getPowerCycle(uint _r) public returns (uint[] memory) {        
        o.push(1);
        o.push(_r);
        while (o[o.length - 1] != 1) {
            o.push((o[o.length - 1] * _r) % MODULUS);
        }
              
        delete o[o.length - 1];
        return o;
    }

    function isPowerOf2(uint _x) internal pure returns (bool) {
        if (_x % 2 == 0) {
            return isPowerOf2(_x.div(2));
        }

        if (_x == 1) {
            return true;
        }

        return false;
    }

    function getPseudorandomIndices(bytes32 _seed, uint _modulus, uint _count, uint _excludeMultiplesOf) internal returns (uint[] memory) {
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

    function _simple_ft(uint[] memory _vals, uint _modulus, uint[] memory _roots) internal returns (uint[] memory) {
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

    function _fft(uint[] memory _vals, uint _modulus, uint[] memory _roots) internal returns (uint[] memory) {
        if (_roots.length <= 4) {
            return _simple_ft(_vals, _modulus, _roots);
        }

        uint halfLength = _vals.length.div(2);

        uint[] memory L = new uint[](halfLength);
        uint[] memory R = new uint[](halfLength);
        uint[] memory _eRoots = new uint[](halfLength);

        uint i;
        for (i = 0; i < _vals.length; i++) {
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

    function fft(uint[] memory _vals, uint _modulus, uint _rootOfUnity, bool _isInv) internal returns (uint[] memory) {
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

        if (_isInv) {
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

    // finite (polynomial) field operations
    function add(uint _x, uint _y) internal returns (uint) {
        return _x.add(_y).mod(MODULUS);
    }

    function sub(uint _x, uint _y) internal returns (uint) {
        return _x.sub(_y).mod(MODULUS);
    }

    function mul(uint _x, uint _y) internal returns (uint) {
        return _x.mul(_y).mod(MODULUS);
    }

    function exp(uint _x, uint _p) internal returns (uint) {
        uint out = 1;
        for (uint i = 1; i < _p; i++) {
            out = out.mul(_x).mod(MODULUS);
        }
        return out;
    }

    function inv(uint _a) internal returns (uint) {
        if (_a == 0) {
            return 0;
        }

        uint lm = 1;
        uint hm = 0;
        uint low = _a.mod(MODULUS);
        uint high = MODULUS;
        while (low > 1) {
            uint r = high.div(low);
            uint nm  = hm.sub(lm.mul(r));
            uint n = high.sub(low.mul(r));

            lm = nm;
            low = n;
            hm = lm;
            high= low;
        }
        return lm.mod(MODULUS);
    }

    function multiInv(uint[] memory _vals) internal returns (uint[] memory) {
        uint[] memory partials = new uint[](_vals.length + 1);
        partials[0] = 1;
        uint i;
        uint y;
        for (i = 0; i < _vals.length; i++) {
            y = 1;
            if (_vals[i] != 0) {
                y = _vals[i];
            }
            partials[i+1] = mul(partials[i], y);
        }

        uint[] memory out = new uint[](_vals.length);
        uint v = partials[partials.length -1];
        for (i = out.length - 1; i >= 0; i--) {
            y = 1;
            if (_vals[i] != 0) {
                out[i] = mul(partials[i], v);
                y = _vals[i];
            } else {
                out[i] = 0;
            }
            v = mul(v, y);
        }
        return out;
    }

    function div(uint _x, uint _y) internal returns (uint) {
        return _x.mul(inv(_y)).mod(MODULUS);
    }

    function evalPolyAt(uint[] memory _p, uint _x) internal returns (uint) {
        uint out = 0;
        uint powerOfX = 1;

        for (uint i = 0; i < _p.length; i++) {
            out = out.add(powerOfX.mul(_p[i]));
            powerOfX = powerOfX.mul(_x).mod(MODULUS);
        }
        return out.mod(MODULUS);
    }

    function addPolys(uint[] memory _a, uint[] memory _b) internal returns (uint[] memory) {
        uint l = _a.length;
        if (l < _b.length) {
            l = _b.length;
        }

        uint[] memory out = new uint[](l);
        for (uint i = 0; i < l; i++) {
            uint aCoff = 0;
            if (i < _a.length) {
                aCoff = _a[i];
            }
            uint bCoff = 0;
            if (i < _b.length) {
                bCoff = _b[i];
            }
            out[i] = aCoff.add(bCoff).mod(MODULUS);
        }
        return out;
    }

    function subPolys(uint[] memory _a, uint[] memory _b) internal returns (uint[] memory) {
        uint l = _a.length;
        if (l < _b.length) {
            l = _b.length;
        }

        uint[] memory out = new uint[](l);
        for (uint i = 0; i < l; i++) {
            uint aCoff = 0;
            if (i < _a.length) {
                aCoff = _a[i];
            }
            uint bCoff = 0;
            if (i < _b.length) {
                bCoff = _b[i];
            }
            out[i] = aCoff.sub(bCoff).mod(MODULUS);
        }
        return out;
    }

    function mulByConst(uint[] memory _a, uint _c) internal returns (uint[] memory) {
        for (uint i = 0; i < _a.length; i++) {
            _a[i] = _a[i].mul(_c).mod(MODULUS);
        }
    }

    function divPolys(uint[] memory _a, uint[] memory _b) internal returns (uint[] memory) {
        require(_a.length >= _b.length);
        if (_b.length == 1) {
            return mulByConst(_a, inv(_b[0]));
        }

        uint apos = _a.length - 1;
        uint bpos = _b.length - 1;
        uint diff = apos - bpos;
        uint[] memory out = new uint[](diff + 1);
        while (diff >= 0) {
            uint quot = div(_a[apos], _b[bpos]);
            out[diff] = quot;
            for (uint i = bpos; i >= 0; i--) {
                _a[diff+i] -= _b[i].mul(quot);
            }
            apos -= 1;
            diff -= 1;
        }

        return out;
    }

    function zPoly(uint[] memory _xs) internal returns (uint[] memory) {
        uint[] memory out = new uint[](_xs.length + 1);
        out[out.length -1] = 1; // top degree

        for (uint i = 0; i < _xs.length; i++) {
            for (uint j = 0; j < i; j++) {
                out[out.length -1 - 1 - j] = out[out.length -1 - j].mul(_xs[i]).mod(MODULUS);
            }
        }
        return out;
    }

    function mulPolys(uint[2] memory _a, uint[2] memory _b) internal returns (uint[] memory) {
        uint[] memory out = new uint[]((_a.length).add(_b.length).sub(1));

        for (uint i = 0; i < _a.length; i++) {
            for (uint j = 0; j < _b.length; j++) {
                out[i+j] = _a[i] * _b[j];
            }
        }
        return out;
    }

    function lagrangeInterp(uint[] memory _xs, uint[] memory _ys) internal returns (uint[] memory) {
        uint[] memory root = zPoly(_xs);
        require(root.length == _ys.length + 1);

        uint[] memory denoms = new uint[](_xs.length);
        uint[] memory tmp = new uint[](2);
        tmp[1] = 1;

        uint i;
        for (i = 0; i < _xs.length; i++) {
            tmp[0] = _xs[i];
            denoms[i] = evalPolyAt(divPolys(root, tmp), _xs[i]);
        }

        uint[] memory invdenoms = multiInv(denoms);
        uint[] memory out = new uint[](_ys.length);

        for (i = 0; i < _xs.length; i++) {
            uint yslice = mul(_ys[i], invdenoms[i]);
            tmp[0] = _xs[i];
            uint[] memory num = divPolys(root, tmp);
            for (uint j = 0; j < _ys.length; j++) {
                if ((num[j] != 0) && (_ys[i] != 0)) {
                    out[j] += num[j].mul(yslice).mod(MODULUS);
                }
            }
        }
        return out;
    }

    // optimized for degree 2
    function lagrangeInterp2(uint[2] memory _xs, uint[2] memory _ys) internal returns (uint[] memory) {
        uint[] memory eq0 = new uint[](2);
        uint[] memory eq1 = new uint[](2);
        eq0[0] = - _xs[1].mod(MODULUS);
        eq1[0] = - _xs[0].mod(MODULUS);
        eq0[1] = 1;
        eq1[1] = 1;
        uint e0 = evalPolyAt(eq0, _xs[0]);
        uint e1 = evalPolyAt(eq1, _xs[1]);
        uint invall = inv(e0.mul(e1));
        uint invY0 = _ys[0].mul(invall).mul(e1);
        uint invY1 = _ys[1].mul(invall).mul(e0);

        uint[] memory out = new uint[](2);
        for (uint i = 0; i < 2; i++) {
            out[i] = ((eq0[i].mul(invY0)).add(eq1[i].mul(invY1))).mod(MODULUS);
        }
        return out;
    }

    // optimized for degree 4
    function lagrangeInterp4(uint[4] memory _xs, uint[4] memory _ys) internal returns (uint[] memory) {
        uint[] memory xxs = new uint[](6);
        xxs[0] = _xs[0].mul(_xs[1]).mod(MODULUS); // 01
        xxs[1] = _xs[0].mul(_xs[2]).mod(MODULUS); // 02
        xxs[2] = _xs[0].mul(_xs[3]).mod(MODULUS); // 03
        xxs[3] = _xs[1].mul(_xs[2]).mod(MODULUS); // 12
        xxs[4] = _xs[1].mul(_xs[3]).mod(MODULUS); // 13
        xxs[5] = _xs[2].mul(_xs[3]).mod(MODULUS); // 23

        uint[4][4] memory eqs;
        eqs[0][0] = -xxs[3].mul(_xs[3]).mod(MODULUS);
        eqs[0][1] = xxs[3].add(xxs[4]).add(xxs[5]);
        eqs[0][2] = -(_xs[1].add(_xs[2].add(_xs[3])));
        eqs[0][3] = 1;

        eqs[1][0] = -xxs[1].mul(_xs[3]).mod(MODULUS);
        eqs[1][1] = xxs[1].add(xxs[2]).add(xxs[5]);
        eqs[1][2] = -(_xs[0].add(_xs[2].add(_xs[3])));
        eqs[1][3] = 1;

        eqs[2][0] = -xxs[0].mul(_xs[3]).mod(MODULUS);
        eqs[2][1] = xxs[0].add(xxs[2]).add(xxs[4]);
        eqs[2][2] = -(_xs[0].add(_xs[1].add(_xs[3])));
        eqs[2][3] = 1;

        eqs[3][0] = -xxs[0].mul(_xs[2]).mod(MODULUS);
        eqs[3][1] = xxs[0].add(xxs[1]).add(xxs[3]);
        eqs[3][2] = -(_xs[0].add(_xs[1].add(_xs[2])));
        eqs[3][3] = 1;

        uint[] memory es = new uint[](4);

        es[0] = evalPolyAt4(eqs[0], _xs[0]);
        es[1] = evalPolyAt4(eqs[1], _xs[1]);
        es[2] = evalPolyAt4(eqs[2], _xs[2]);
        es[3] = evalPolyAt4(eqs[3], _xs[3]);

        uint i = inv(es[0].mul(es[1]).mul(es[2].mul(es[3])));

        uint[] memory invYs = new uint[](4);
        invYs[0] = _ys[0].mul(i).mul(es[1]).mul(es[2].mul(es[3])).mod(MODULUS);
        invYs[1] = _ys[1].mul(i).mul(es[0]).mul(es[2].mul(es[3])).mod(MODULUS);
        invYs[2] = _ys[2].mul(i).mul(es[3]).mul(es[0].mul(es[1])).mod(MODULUS);
        invYs[3] = _ys[3].mul(i).mul(es[2]).mul(es[0].mul(es[1])).mod(MODULUS);

        uint[] memory out = new uint[](4);
        for (i = 0; i < 4; i++) {
            out[i] = ((eqs[0][i].mul(invYs[0])
            ).add(eqs[1][i].mul(invYs[1])
            ).add(eqs[2][i].mul(invYs[2])
            ).add(eqs[3][i].mul(invYs[3])
            )).mod(MODULUS);
        }
        return out;
    }

    function multiInterp4(uint[4][] memory _xsets, uint[4][] memory _ysets) internal returns (uint[][] memory) {
        require(_xsets.length == _ysets.length);
        uint[][] memory out = new uint[][](_xsets.length);

        for (uint i = 0; i < _xsets.length; i++) {
            out[i] = lagrangeInterp4(_xsets[i], _ysets[i]);
        }
        return out;
    }

    function evalPolyAt4(uint[4] memory _p, uint _x) internal returns (uint) {
        uint out = 0;
        uint powerOfX = 1;

        for (uint i = 0; i < _p.length; i++) {
            out = out.add(powerOfX.mul(_p[i]));
            powerOfX = powerOfX.mul(_x).mod(MODULUS);
        }
        return out.mod(MODULUS);
    }

    function evalQuartic(uint[] memory p, uint x) internal returns (uint) {
        require(p.length == 4);
        uint xsq = x.mul(x).mod(MODULUS);
        uint xcb = xsq.mul(x).mod(MODULUS);
        return (p[0] + p[1].mul(x) + p[2] * xsq + p[3] * xcb).mod(MODULUS);
    }
}
