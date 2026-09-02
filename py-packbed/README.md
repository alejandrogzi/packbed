# py-packbed

Python bindings for `packbed`, a Rust library that groups BED records into
overlapping transcript components.

```python
import packbed

components = packbed.pack(
    ["reference.bed", "query.bed"],
    ["reference", "query"],
    overlap_type="exon",
)
```

The Python package intentionally exposes only conversion of the Rust packed
component result into Python objects. It does not provide BED file writing,
serialization, colorization, or CLI behavior.
