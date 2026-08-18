# Expected metrics

No accepted numerical baseline has been established for this legacy case yet.

After a successful run has been reviewed, add `metrics.json` containing reference values and tolerances. Do not derive an accepted baseline automatically from a failed or unreviewed run.

Example shape:

```json
{
  "elastic_E": {
    "x": { "value": 0.0, "atol": 1e-8, "rtol": 1e-6 },
    "fun": { "value": 0.0, "atol": 1e-10, "rtol": 1e-6 }
  }
}
```
