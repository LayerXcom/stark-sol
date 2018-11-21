pragma solidity ^0.4.24;

import "./SafeMath.sol";
import "./BytesLib.sol";

library Merkle {
    using SafeMath for uint;
    using BytesLib for bytes;

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
        bytes[] _proof
    ) 
        internal 
        pure 
        returns (bytes)
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

    function merkelize(uint[]) internal returns (bytes32[]) {
        bytes32[] memory a = new bytes32[](3);
        return a;
    }

    function permute4(uint[]) internal returns (uint[]) {
        uint[] memory a = new uint[](3);
        return a;
    }
        
}