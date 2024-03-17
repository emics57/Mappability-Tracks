```mappabilityTracks.py```: manually finds the mappable regions by iterating through all read coordinates

```mappabilityBEDtools.py```: outputs read coordinates as a BED file and uses bedtools merge to merge reads together to get mappable regions

Both of these scripts (should) do the same thing and output the same exact regions in the file. It's just that the latter requires BEDtools. The bedtools version was created as a 'back up/sanity check' for mappabilityTracks.py
