# Supported reference modes

## Three-reference plant mode
Nuclear + chloroplast + mitochondrial references. The workflow can report inclusive organelle-compatible, strict organelle-compatible, and nuclear–organelle multi-reference fractions.

## Two-reference mode
Nuclear + one organelle reference (used for the chicken validation). The same set logic is applied to the available organelle.

## Organelle-only mode
If no nuclear reference is available, the workflow can report inclusive compatibility to the supplied organelle reference(s). It intentionally does not fabricate `strict` or `multi-reference` values, because those quantities require a nuclear comparison.

This mode is useful for non-model organisms or archive audits in which only organelle references are available.
