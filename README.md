# ChemCards - medicinal chemistry structure flashcards

A Streamlit flashcard app for drilling chemical structures of drugs,
neurotransmitters, amino acids, and nucleic acid bases. Structures render from
SMILES with RDKit; formula and MW are derived automatically.

## Run locally

```bash
pip install -r requirements.txt
streamlit run app.py
```

## Decks: daily regimen and free practice

The sidebar has two collapsible decks plus shared options.

**My daily deck.** Pick the categories to drill and a number of cards per day,
then "Start today's session" — the deck is a fresh random sample of that size.
Work through every card (answer or give up) and a "Finish today's session"
button appears; finishing shows a celebration splash and advances your streak.
Category choice, cards/day, and the streak are remembered between visits via a
browser cookie (`streamlit-cookies-controller`), so you can return the next day
and continue the streak. The streak increments on consecutive days, doesn't
double-count two finishes in one day, and resets after a missed day. Settings
are editable any time in the same expander.

**Deck.** Free practice — the full category toggles (with All Biochemistry /
All Drugs master switches), shuffle, and "Review missed only"; click "Train this
deck" to switch to it. You can free-practice before, after, or instead of the
daily session; starting/resuming the daily deck never disturbs the free
selection and vice versa.

Shared below both: study direction and the formula/MW hint toggle.

## Study modes

**Structure -> Name.** The structure is shown; you type the name and click
Check. Matching is case/space/hyphen-insensitive and accepts aliases drawn from
the `abbrev` field (three-letter and one-letter codes, brand/INN names, e.g.
`Phe`, `F`, `adrenaline`, `Prozac`). A near-miss prompts you to fix spelling
rather than failing you outright. You keep trying until correct or you reveal.

**Name -> Structure.** The name is shown; you draw the molecule in an embedded
Ketcher editor and click Apply. Your drawing is scored by:

- **Empirical formula table** - per-element counts (incl. explicit H via
  `AddHs`) for reference vs. your guess, with the delta; rows are green where
  they match, red where they don't.
- **Match level** - *exact* (constitution + stereochemistry), *constitution*
  (graph matches, stereo differs - counted as correct), or a miss.
- **Descriptor deltas** - heavy atoms, H-bond donors/acceptors, aromatic and
  aliphatic rings, phrased as "N too many / N missing".
- **ECFP4 Tanimoto** to the target.
- **MCS highlight** on your drawing - two passes after the Moleculardle style:
  a tight pass (`matchValences`, `ringMatchesRingOnly`) highlighted blue, and a
  loose pass (`bondCompare=CompareAny`) highlighted red where the tight pass
  missed. In rings, a blue atom/bond drops to red if its aromaticity disagrees
  with the target. (Note: the upstream `IsInRing == True` comparison was a bug -
  it compared a bound method to True and never fired; fixed here to `IsInRing()`.)

Both modes share the sidebar: toggle compound groups, shuffle, session score,
and a "Review missed only" deck built from cards you missed this session.

## Stereochemistry

Drawing exact stereochemistry in Ketcher is tedious, so a correct 2D
constitution counts as solved; the app tells you when stereo differs from the
reference. If you want to *require* stereo, change the `match_level` acceptance
in `app.py` to treat only `"exact"` as solved.

## Adding compounds

Edit `compounds.json`. Each entry uses a `categories` **list**, so a compound
can belong to several groups at once:

```json
{ "name": "Tramadol", "smiles": "CN(C)CC1CCCCC1(O)c1cccc(OC)c1",
  "categories": ["Drugs", "Opioids", "SSRIs/SNRIs"], "abbrev": "Ultram" }
```

- A compound appears in the deck if **any** of its categories is enabled; it's
  shown once regardless of how many active groups it belongs to.
- A new category string automatically becomes a new sidebar toggle. The display
  order is set by `_ORDER` in `app.py`; unlisted categories sort in after it.
- The legacy single `"category": "..."` form is still accepted and treated as a
  one-element list.
- Formula/MW are computed at load; invalid SMILES are skipped.
- The `abbrev` field doubles as the alias list for name-mode answers (split on
  `/` and `,`), so it's where brand names and alternate names go (e.g. glutamate
  is accepted for glutamic acid; Prozac/Sarafem for fluoxetine).

The deck ships as a validated build: `build.py` pairs every drug SMILES with an
expected molecular formula and refuses to emit anything that fails to parse or
whose RDKit-computed formula disagrees. Re-run `python build.py` after edits to
regenerate `compounds.json` with the same check. Structures are neutral
free-base/free-acid parent forms (not salts).

### Categories in the shipped deck

Biomolecules: Amino acids (20), Neurotransmitters (11), Nucleic acid bases (5).
Glycine, aspartic acid, and glutamic acid are tagged as both amino acids and
neurotransmitters. Drug classes (all also tagged `Drugs`, the 149-compound
catch-all): Antipsychotics (16), SSRIs/SNRIs (14), Sedatives & anxiolytics (20),
Opioids (15), Stimulants (10), NSAIDs & analgesics (12), Statins (7),
Antihistamines & allergy (16), Antibiotics (16), Antifungals (11), Others (13,
for drugs that don't fit a class above — PPIs, anticoagulants, antihypertensives,
a bronchodilator, an antiviral, etc.). Agrochemicals, which are **not** tagged
`Drugs`: Insecticides (18), Herbicides (18). 218 compounds total. Overlaps are
allowed where appropriate (e.g. tramadol is an opioid and an SNRI).

## Deploy on Streamlit Community Cloud

Push the folder to GitHub and point an app at `app.py`; `requirements.txt`
installs RDKit, Ketcher, and pandas automatically.
