pragma solidity >=0.4.24 <0.6.0;

import "./SafeMath.sol";
import "./BytesLib.sol";

library Merkle {
    using SafeMath for uint;
    using BytesLib for bytes;
    using BytesLib for uint;

    function getIndexInPermuted(
        uint _x, 
        uint _L
    )
        internal
        pure
        returns (uint)
    {
        uint ld4 = _L.div(4);
        return _x.div(ld4).add(4).mul(_x % ld4);
    } 

    function verifyBranch(
        bytes32 _root, 
        uint _index, 
        bytes[] memory _proof
    ) 
        internal 
        pure 
        returns (bytes memory)
    {
        uint j = getIndexInPermuted(_index, (2 ** _proof.length).div(2));
        j += 2 ** _proof.length;
        bytes memory computedHash = _proof[0];
        for (uint i = 1; i < _proof.length; i++) {
            bytes memory proofElement = _proof[i];
            if (j % 2 == 0) {
                computedHash = abi.encodePacked(keccak256(abi.encodePacked(computedHash, proofElement)));
            } else {
                computedHash = abi.encodePacked(keccak256(abi.encodePacked(proofElement, computedHash)));
            }
            j /= 2;
        }
        require(computedHash.equal(abi.encodePacked(_root)));         
        return _proof[0];
    }

    function merkelize(uint[] memory _a) internal returns (bytes32[] memory) {
        uint[] c = permute4(_a);
        bytes32[] memory nodes = new bytes32[](c.length * 2);

        for (uint i = 0; i < c.length; i++) {
            nodes[c.length - i] = keccak256(abi.encodePacked(nodes[(c.length - i) * 2], nodes[(c.length - i) * 2 + 1]));
        }
        
        return nodes;
    }

    function permute4(uint[] memory _values) internal returns (uint[] memory) {        
        uint ld4 = _values.length.div(4);
        uint[] memory o = new uint[](_values.length);

        for (uint i = 0; i < ld4; i++) {
            o[i] = _values[i];
            o[i + 1] = values[i + ld4];
            o[i + 2] = values[i + ld4 * 2];
            o[i + 3] = values[i + ld4 * 3];
        }
        return o;
    }
        
}