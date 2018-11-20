# stark-sol
a solidity implementation of verification in a STARK on a MIMC calculation

## setup

```
git clone git@github.com:LayerXcom/stark-sol.git
cd stark-sol
pip install -r requirements.txt
yarn install
```


## test
generate a proof for stark. It'll output as a `.json` file, `stark_proof.json`.
```
python -m pytest test.py
```

truffle test for the verifier contract
```
truffle development
> compile
> migrate
> test
```