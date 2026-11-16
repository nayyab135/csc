# Notes for Claude

- Always reply to the user in English, even when the user writes in Hindi or Urdu. This applies to chat replies and to every document produced (PDFs, notes, code comments). Use simple, plain wording.
- The user runs the MATLAB codes on their own machine (MATLAB with parpool); there is no MATLAB in the cloud session. Octave can be installed to parse-check and smoke-test, but local functions must be moved to the top of a test copy for Octave.
- When the user asks for figures only, give a plotting-only script that uses the workspace variables; never `clear` the workspace in such a script.
