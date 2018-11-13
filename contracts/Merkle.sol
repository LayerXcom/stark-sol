pragma solidity ^0.4.24;

import "./SafeMath.sol";

library Merkle {
    using SafeMath for uint;

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
        bytes32[] _proof
    ) 
        internal 
        pure 
        returns (bytes32)
    {
        uint j = getIndexInPermuted(_index, (2 ** _proof.length).div(2));
        j += 2 ** _proof.length;
        bytes32 computedHash = _proof[0];
        for (uint i = 1; i < _proof.length; i++) {
            bytes32 proofElement = _proof[i];
            if (j % 2 == 0) {
                computedHash = keccak256(abi.encodePacked(computedHash, proofElement));
            } else {
                computedHash = keccak256(abi.encodePacked(proofElement, computedHash));
            }
            j /= 2;
        }
        require(computedHash == _root);
        return _proof[0];
    }
}