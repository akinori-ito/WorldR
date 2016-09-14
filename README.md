# WorldR: An R wrapper of World 

## Overview
WorldR is a wrapper of World, the speech analysis and synthesis system.
- [World development page on github](https://github.com/mmorise/World)

It provides the following functionalities:
- Analyze a speech signal (the Wave object) into F0, spectrum and aperiodicity
- Synthesize a speech signal from the F0, spectrum and aperiodicity
- Change F0 of the speech without changing the tempo
- Stretch and compress the speech signal without changing F0
- Change the vocal identity (virtually the length of the vocal tract) without changing F0 and tempo

## Dependency
- tuneR
