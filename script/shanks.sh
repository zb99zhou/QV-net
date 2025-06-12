#!/bin/bash

echo "[Table 6]"
cd ..
cargo test test_shanks --release -- --nocapture 2>&1 | grep -E "Time elapsed in Shanks for scalar 2\^[0-9]+: [0-9.]+ms"
