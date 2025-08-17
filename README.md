# OFDM Projects

This repository contains MATLAB projects focused on **Orthogonal Frequency Division Multiplexing (OFDM)** and its enhancements.  
All simulations use **16-QAM modulation** and include numerical outputs as well as plots.

## 📂 Contents

### 1. Basic OFDM
- **File:** `ofdm.m`  
- **Description:**  
  Implements a baseline OFDM system with 16-QAM modulation over an AWGN channel.  
  Useful as a starting point for comparison with more advanced schemes.

---

### 2. OFDM with LDPC
- **Files:**  
  - `ofdm_ldpc_random_number.m`  
  - `ofdm_ldpc_foto.m`  
- **Description:**  
  Adds **LDPC channel coding** to the OFDM system.  
  Two variants are provided:  
  - `random_number` → Input bits generated randomly.  
  - `foto` → Transmission of an image file encoded with OFDM+LDPC.  

---

### 3. OFDM with Multipath & Doppler Fading
- **Files:**  
  - `ofdm_ldpc_pilotbased_mudofading_16qam_TOP.m`  
  - `ofdm_pilotbased_mudofading_16qam_TOP.m`  
- **Description:**  
  Models an **OFDM system over a multipath fading channel with Doppler effect**.  
  - Uses **pilot-based channel estimation**.  
  - Includes **LDPC coding** for error correction.  
  - Evaluates robustness under time-varying channel conditions.  

---

## 🛠 Requirements
- MATLAB R2024b (tested)  
- Communications System Toolbox (for LDPC and channel models)  

---

## 📑 Notes
These codes are part of my thesis work on **multicarrier modulation techniques** for digital communication systems.  
Focus is placed on:  
- Comparing uncoded vs LDPC-coded OFDM.  
- Studying the effects of multipath and Doppler fading.  
- Evaluating pilot-based channel estimation.
