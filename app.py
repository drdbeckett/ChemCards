"""
ChemCards - medicinal chemistry structure flashcards.
Run with:  streamlit run app.py

Two interaction modes:
  - Structure -> Name : structure is shown, you TYPE the name (graded).
  - Name -> Structure : name is shown, you DRAW it in Ketcher; the app scores
                        your drawing by empirical formula, descriptor deltas,
                        Tanimoto, and a blue/red maximum-common-substructure
                        highlight (in the Moleculardle drawing style).
"""

import json
import random
import difflib
from pathlib import Path
from datetime import date, timedelta
from collections import Counter, defaultdict

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
from streamlit_cookies_controller import CookieController
from streamlit_ketcher import st_ketcher
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, rdMolDescriptors, rdFMCS, rdFingerprintGenerator
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor

rdDepictor.SetPreferCoordGen(True)   # cleaner CoordGen 2D depictions

DATA_PATH = Path(__file__).parent / "compounds.json"
st.set_page_config(page_title="ChemCards", page_icon="\U0001F9EA", layout="wide")

# transparent blue (tight match) / red (loose match), Moleculardle palette
_BLUE = (0.0, 0.0, 1.0, 0.25)
_RED = (1.0, 0.0, 0.0, 0.25)
_FPGEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


# ---------------------------------------------------------------------------
# Data + chemistry helpers (cached)
# ---------------------------------------------------------------------------
@st.cache_data
def load_compounds():
    raw = json.loads(DATA_PATH.read_text())
    out = []
    for c in raw:
        mol = Chem.MolFromSmiles(c["smiles"])
        if mol is None:
            continue
        c = dict(c)
        # Accept either a "categories" list or a legacy single "category" string.
        if "categories" not in c:
            c["categories"] = [c["category"]] if c.get("category") else []
        c["canonical"] = Chem.MolToSmiles(mol)
        c["formula"] = rdMolDescriptors.CalcMolFormula(mol)
        c["mw"] = round(Descriptors.MolWt(mol), 2)
        c["aliases"] = _aliases(c)
        out.append(c)
    return out


def _norm(s: str) -> str:
    return "".join(s.lower().split()).replace("-", "")


def _aliases(card) -> set:
    """Acceptable answers: the name plus any token in `abbrev` (split on / and ,)."""
    al = {_norm(card["name"])}
    for tok in card.get("abbrev", "").replace(",", "/").split("/"):
        tok = tok.strip()
        if tok:
            al.add(_norm(tok))
    return al


@st.cache_data
def render_structure(smiles: str, size: int = 420) -> bytes:
    mol = Chem.MolFromSmiles(smiles)
    d = rdMolDraw2D.MolDraw2DCairo(size, size)
    d.drawOptions().bondLineWidth = 2
    rdMolDraw2D.PrepareAndDrawMolecule(d, mol)
    d.FinishDrawing()
    return d.GetDrawingText()


@st.cache_data
def render_comparison(target_smiles: str, guess_smiles: str, size: int = 420) -> bytes:
    """Highlight the guess by its maximum common substructure with the target.

    Two MCS passes (after Moleculardle):
      tight  - matchValences + ringMatchesRingOnly  -> blue
      loose  - bondCompare=CompareAny               -> red (only where tight missed)
    In rings, a blue atom/bond is downgraded to red if its aromaticity disagrees
    with the target.
    """
    tm = Chem.MolFromSmiles(target_smiles)
    gm = Chem.MolFromSmiles(guess_smiles)

    loose = rdFMCS.FindMCS([tm, gm], bondCompare=rdFMCS.BondCompare.CompareAny)
    lmol = Chem.MolFromSmarts(loose.smartsString)
    lg = gm.GetSubstructMatch(lmol)
    tight = rdFMCS.FindMCS([tm, gm], matchValences=True, ringMatchesRingOnly=True)
    tmol = Chem.MolFromSmarts(tight.smartsString)
    tg = gm.GetSubstructMatch(tmol)
    tt = tm.GetSubstructMatch(tmol)

    ath, arad = defaultdict(list), {}
    for sid in range(len(tg)):
        gid, tid = tg[sid], tt[sid]
        arad[gid] = 0.3
        if gm.GetAtomWithIdx(gid).IsInRing():
            same_arom = (gm.GetAtomWithIdx(gid).GetIsAromatic()
                         == tm.GetAtomWithIdx(tid).GetIsAromatic())
            ath[gid].append(_BLUE if same_arom else _RED)
        else:
            ath[gid].append(_BLUE)
    for gid in lg:
        if gid not in ath:
            ath[gid].append(_RED)

    bnd = defaultdict(list)
    for b in tmol.GetBonds():
        g1, g2 = tg[b.GetBeginAtomIdx()], tg[b.GetEndAtomIdx()]
        t1, t2 = tt[b.GetBeginAtomIdx()], tt[b.GetEndAtomIdx()]
        gb = gm.GetBondBetweenAtoms(g1, g2)
        tb = tm.GetBondBetweenAtoms(t1, t2)
        if gb.IsInRing():
            same_arom = gb.GetIsAromatic() == tb.GetIsAromatic()
            bnd[gb.GetIdx()].append(_BLUE if same_arom else _RED)
        else:
            bnd[gb.GetIdx()].append(_BLUE)
    for b in lmol.GetBonds():
        g1, g2 = lg[b.GetBeginAtomIdx()], lg[b.GetEndAtomIdx()]
        gb = gm.GetBondBetweenAtoms(g1, g2).GetIdx()
        if gb not in bnd:
            bnd[gb].append(_RED)

    d = rdMolDraw2D.MolDraw2DCairo(size, size)
    d.drawOptions().bondLineWidth = 2
    d.DrawMoleculeWithHighlights(gm, "", dict(ath), dict(bnd), arad, {})
    d.FinishDrawing()
    return d.GetDrawingText()


def element_counts(mol) -> Counter:
    return Counter(a.GetSymbol() for a in Chem.AddHs(mol).GetAtoms())


def formula_table(target_mol, guess_mol):
    tc, gc = element_counts(target_mol), element_counts(guess_mol)
    elems = sorted(set(tc) | set(gc), key=lambda e: (e != "C", e != "H", e))
    df = pd.DataFrame([
        {"Element": e, "Reference": tc.get(e, 0), "Yours": gc.get(e, 0),
         "\u0394 (yours\u2212ref)": gc.get(e, 0) - tc.get(e, 0)}
        for e in elems
    ])

    def _row_style(row):
        ok = row["\u0394 (yours\u2212ref)"] == 0
        bg = "rgba(40,170,80,0.18)" if ok else "rgba(210,60,60,0.18)"
        return [f"background-color: {bg}"] * len(row)

    return df.style.apply(_row_style, axis=1)


def tanimoto(a_mol, b_mol) -> float:
    return DataStructs.TanimotoSimilarity(
        _FPGEN.GetFingerprint(a_mol), _FPGEN.GetFingerprint(b_mol))


def match_level(target_smiles: str, guess_mol) -> str:
    """'exact' (incl. stereo), 'constitution' (graph matches, stereo differs), or 'no'."""
    tm = Chem.MolFromSmiles(target_smiles)
    if Chem.MolToSmiles(guess_mol) == Chem.MolToSmiles(tm):
        return "exact"
    g2 = Chem.MolFromSmiles(Chem.MolToSmiles(guess_mol))
    t2 = Chem.MolFromSmiles(target_smiles)
    Chem.RemoveStereochemistry(g2)
    Chem.RemoveStereochemistry(t2)
    return "constitution" if Chem.MolToSmiles(g2) == Chem.MolToSmiles(t2) else "no"


def property_readout(num: int, prop: str) -> str:
    num = int(num)
    if num == 0:
        return f"\u2713 correct number of {prop}s"
    if num > 0:
        return f"{num} {prop}{'s' if num != 1 else ''} too many"
    return f"{-num} {prop}{'s' if num != -1 else ''} missing"


def descriptor_deltas(target_mol, guess_mol) -> list:
    pairs = [
        (rdMolDescriptors.CalcNumHeavyAtoms, "heavy atom"),
        (rdMolDescriptors.CalcNumHBD, "H-bond donor"),
        (rdMolDescriptors.CalcNumHBA, "H-bond acceptor"),
        (rdMolDescriptors.CalcNumAromaticRings, "aromatic ring"),
        (rdMolDescriptors.CalcNumAliphaticRings, "aliphatic ring"),
    ]
    return [property_readout(fn(guess_mol) - fn(target_mol), label) for fn, label in pairs]


COMPOUNDS = load_compounds()
_all_cats = {cat for c in COMPOUNDS for cat in c["categories"]}

# Sidebar groupings. "Drugs" is the meta-tag every drug also carries; it isn't a
# filter checkbox of its own anymore (the subclasses below cover every drug), so
# it's omitted here. "Others" is the drug-class catch-all, shown as "Other Drugs".
BIO_CATS = [c for c in ["Amino acids", "Neurotransmitters", "Nucleobases"]
            if c in _all_cats]
DRUG_CATS = [c for c in ["Antipsychotics", "SSRIs/SNRIs", "Sedatives & anxiolytics",
                         "Opioids", "Stimulants", "NSAIDs & analgesics", "Statins",
                         "Antihistamines & allergy", "Antibiotics", "Antifungals",
                         "Others"] if c in _all_cats]
AGRO_CATS = [c for c in ["Insecticides", "Herbicides"] if c in _all_cats]
CATEGORIES = BIO_CATS + DRUG_CATS + AGRO_CATS
DISPLAY = {"Others": "Other Drugs"}


# ---------------------------------------------------------------------------
# Session state
# ---------------------------------------------------------------------------
def init_state():
    ss = st.session_state
    ss.setdefault("deck", [])
    ss.setdefault("pos", 0)
    ss.setdefault("deck_key", None)
    ss.setdefault("nonce", 0)        # bumps every card change -> fresh widgets
    ss.setdefault("seen", 0)
    ss.setdefault("correct", 0)
    ss.setdefault("missed", set())
    # per-card scratch
    ss.setdefault("answered", False)
    ss.setdefault("solved", False)
    ss.setdefault("attempts", 0)
    ss.setdefault("last_guess", None)
    ss.setdefault("last_processed", None)
    # daily regimen / streak
    ss.setdefault("active_mode", "free")     # "free" or "daily"
    ss.setdefault("daily_answered", set())   # compound indices answered today
    ss.setdefault("daily_completed_flag", False)
    ss.setdefault("show_splash", None)       # streak number to celebrate, or None
    ss.setdefault("finish_pending", False)   # finish requested; write cookie next run


def reset_card():
    ss = st.session_state
    ss.answered = False
    ss.solved = False
    ss.attempts = 0
    ss.last_guess = None
    ss.last_processed = None
    ss.nonce += 1


def build_deck(active_cats, shuffle, missed_only):
    ss = st.session_state
    pool = [i for i, c in enumerate(COMPOUNDS)
            if any(cat in active_cats for cat in c["categories"])
            and (not missed_only or i in ss.missed)]
    if shuffle:
        random.shuffle(pool)
    ss.deck = pool
    ss.pos = 0
    reset_card()


def advance(step):
    ss = st.session_state
    if ss.deck:
        ss.pos = (ss.pos + step) % len(ss.deck)
        reset_card()


def grade(solved: bool):
    """Count the current card exactly once."""
    ss = st.session_state
    if ss.answered:
        return
    ss.answered = True
    ss.solved = solved
    ss.seen += 1
    idx = ss.deck[ss.pos]
    if solved:
        ss.correct += 1
        ss.missed.discard(idx)
    else:
        ss.missed.add(idx)
    if ss.active_mode == "daily":
        ss.daily_answered.add(idx)


# ---------------------------------------------------------------------------
# Daily regimen + streak (persisted in a cookie)
# ---------------------------------------------------------------------------
COOKIE_KEY = "chemcards"
DEFAULT_SETTINGS = {"daily_cats": [], "daily_n": 20, "streak": 0,
                    "last_completed": None}
# Show the Storage debug panel only when enabled. Set True here, or append
# ?debug=1 to the app URL to toggle it on without editing the file.
DEBUG = False


def update_streak(s):
    """Advance the streak on completion. Same day = no change; consecutive day
    = +1; any gap = reset to 1. Returns the new streak."""
    today = date.today().isoformat()
    yesterday = (date.today() - timedelta(days=1)).isoformat()
    last = s.get("last_completed")
    if last == today:
        pass
    elif last == yesterday:
        s["streak"] = int(s.get("streak", 0)) + 1
        s["last_completed"] = today
    else:
        s["streak"] = 1
        s["last_completed"] = today
    return int(s["streak"])


def active_streak(s):
    """The streak only counts as 'alive' if the last completion was today or
    yesterday; otherwise it's broken (shown as 0 until a new session)."""
    last = s.get("last_completed")
    alive = last in (date.today().isoformat(),
                     (date.today() - timedelta(days=1)).isoformat())
    return int(s.get("streak", 0)) if alive else 0


def daily_done_today(s):
    return s.get("last_completed") == date.today().isoformat()


def start_daily():
    ss = st.session_state
    cats = ss.settings.get("daily_cats", [])
    n = max(1, int(ss.settings.get("daily_n", 20)))
    pool = [i for i, c in enumerate(COMPOUNDS)
            if any(cat in cats for cat in c["categories"])]
    random.shuffle(pool)
    ss.deck = pool[:n]
    ss.pos = 0
    ss.active_mode = "daily"
    ss.daily_answered = set()
    ss.daily_completed_flag = False
    ss.deck_key = None        # force a clean rebuild if we later return to free
    reset_card()


def finish_daily():
    ss = st.session_state
    new = update_streak(ss.settings)
    save_settings()
    ss.show_splash = new
    ss.daily_completed_flag = True


def render_splash(streak):
    ss = st.session_state
    st.balloons()
    msg = ("You've started a new streak!" if streak <= 1
           else f"You have a {streak} day streak!")
    st.markdown(
        "<div style='text-align:center;padding:3rem 1rem'>"
        "<div style='font-size:5rem;line-height:1'>\U0001F525</div>"
        "<h1 style='margin:0.5rem 0'>Congratulations!</h1>"
        f"<h2 style='margin:0;font-weight:400'>{msg}</h2></div>",
        unsafe_allow_html=True,
    )
    if st.button("Continue to free practice \u2192", type="primary",
                 use_container_width=True):
        ss.show_splash = None
        ss.active_mode = "free"
        ss.deck_key = None
        st.rerun()


init_state()
ss = st.session_state

# Cookie-backed settings (daily categories, cards/day, streak, last completion).
controller = CookieController()


def load_settings():
    raw = controller.get(COOKIE_KEY)
    ss.setdefault("settings", dict(DEFAULT_SETTINGS))
    if not ss.get("_cookie_synced") and raw:
        try:
            data = json.loads(raw) if isinstance(raw, str) else raw
            if isinstance(data, dict):
                ss.settings = {**DEFAULT_SETTINGS, **data}
        except Exception:
            pass
        ss._cookie_synced = True


def save_settings():
    try:
        # Default cookie expiry in this library is 1 day, which would drop the
        # streak between sessions; keep it for a year instead.
        controller.set(COOKIE_KEY, json.dumps(ss.settings),
                       max_age=60 * 60 * 24 * 365)
    except Exception:
        pass


load_settings()


# ---------------------------------------------------------------------------
# Sidebar
# ---------------------------------------------------------------------------
def _count(cat):
    return sum(cat in c["categories"] for c in COMPOUNDS)


def _group_count(cats):
    cset = set(cats)
    return sum(bool(cset.intersection(c["categories"])) for c in COMPOUNDS)


with st.sidebar:
    st.title("\U0001F9EA ChemCards")

    # ---------------- My daily deck ----------------
    streak_now = active_streak(ss.settings)
    done_today = daily_done_today(ss.settings)
    daily_label = "\U0001F525 My daily deck"
    if streak_now:
        daily_label += f"  \u2014  {streak_now}\U0001F525"
    with st.expander(daily_label, expanded=not done_today):
        if streak_now:
            st.caption(f"Current streak: {streak_now} "
                       f"day{'s' if streak_now != 1 else ''}"
                       + ("  \u00b7  done today \u2705" if done_today else ""))
        else:
            st.caption("No active streak \u2014 finish a session to start one.")

        # Auto-populate last session's choices. The cookie loads asynchronously,
        # so we can't rely on the widget `default=` (a keyed widget locks in its
        # first value before the cookie arrives). Instead, push the saved values
        # into the widgets' state once, the moment the cookie has synced.
        st.session_state.setdefault("daily_cats_widget", [])
        st.session_state.setdefault("daily_n_widget", 20)
        if ss.get("_cookie_synced") and not ss.get("_daily_hydrated"):
            st.session_state["daily_cats_widget"] = [
                c for c in ss.settings.get("daily_cats", []) if c in CATEGORIES]
            st.session_state["daily_n_widget"] = int(ss.settings.get("daily_n", 20))
            ss._daily_hydrated = True

        sel = st.multiselect("Categories to drill", options=CATEGORIES,
                             format_func=lambda c: DISPLAY.get(c, c),
                             key="daily_cats_widget")
        n_day = st.number_input("Cards per day", min_value=5, max_value=300,
                                step=5, key="daily_n_widget")
        if sel != ss.settings.get("daily_cats") or int(n_day) != int(ss.settings.get("daily_n", 20)):
            ss.settings["daily_cats"] = sel
            ss.settings["daily_n"] = int(n_day)
            save_settings()

        avail = sum(any(cat in sel for cat in c["categories"]) for c in COMPOUNDS)
        start_label = ("Resume today's session" if ss.active_mode == "daily"
                       and not ss.daily_completed_flag else "Start today's session")
        if st.button(start_label, type="primary", use_container_width=True,
                     disabled=not sel):
            start_daily()
            st.rerun()
        st.caption(f"{min(int(n_day), avail)} cards from {avail} available."
                   if sel else "Pick at least one category.")

    # ---------------- Free deck ----------------
    with st.expander("Deck", expanded=ss.active_mode == "free"):
        for _c in CATEGORIES:
            st.session_state.setdefault(f"cat_{_c}", True)
        st.session_state.setdefault("all_bio", True)
        st.session_state.setdefault("all_drugs", True)

        def _apply_master(master_key, cats):
            on = st.session_state[master_key]
            for c in cats:
                st.session_state[f"cat_{c}"] = on

        def _cat_box(cat):
            st.checkbox(f"{DISPLAY.get(cat, cat)} ({_count(cat)})", key=f"cat_{cat}")

        st.caption("Biochemistry")
        st.checkbox(f"**All Biochemistry** ({_group_count(BIO_CATS)})", key="all_bio",
                    on_change=_apply_master, args=("all_bio", BIO_CATS))
        for cat in BIO_CATS:
            _cat_box(cat)

        st.caption("Drugs")
        st.checkbox(f"**All Drugs** ({_group_count(DRUG_CATS)})", key="all_drugs",
                    on_change=_apply_master, args=("all_drugs", DRUG_CATS))
        for cat in DRUG_CATS:
            _cat_box(cat)

        st.caption("Others")
        for cat in AGRO_CATS:
            _cat_box(cat)

        active_cats = [c for c in CATEGORIES if st.session_state[f"cat_{c}"]]
        st.divider()
        shuffle = st.checkbox("Shuffle", value=True)
        missed_only = st.checkbox("Review missed only", value=False,
                                  disabled=not ss.missed)
        if st.button("Train this deck", use_container_width=True):
            ss.active_mode = "free"
            ss.deck_key = None
            st.rerun()

    # ---------------- Global options ----------------
    st.divider()
    mode = st.radio("Study direction",
                    ["Structure \u2192 Name", "Name \u2192 Structure"],
                    help="Structure->Name: type the name. "
                         "Name->Structure: draw it in Ketcher.")
    show_hints = st.checkbox("Show formula + MW with answer", value=True)
    st.metric("Session score", f"{ss.correct}/{ss.seen}" if ss.seen else "0/0")
    if st.button("Clear stats", use_container_width=True):
        ss.seen = ss.correct = 0
        ss.missed = set()

    # ---------------- Storage debug (hidden unless enabled) ----------------
    if DEBUG or st.query_params.get("debug") == "1":
        with st.expander("Storage debug", expanded=False):
            st.write("cookie synced:", bool(ss.get("_cookie_synced")))
            try:
                st.write("cookies seen by app:", controller.getAll())
            except Exception as e:  # noqa: BLE001
                st.write("getAll error:", repr(e))
            st.write("settings in use:", ss.settings)
            c1, c2 = st.columns(2)
            if c1.button("Write test value", use_container_width=True):
                ss.settings["_probe"] = date.today().isoformat()
                save_settings()
                st.toast("Wrote a test value to the cookie.")
            if c2.button("Re-read cookies", use_container_width=True):
                try:
                    controller.refresh()
                except Exception as e:  # noqa: BLE001
                    st.write("refresh error:", repr(e))
                ss._cookie_synced = False
                st.rerun()

# Auto-rebuild the free deck from the toggles only while in free mode, so a
# daily session in progress is never clobbered.
if ss.active_mode == "free":
    deck_key = (tuple(active_cats), shuffle, missed_only)
    if ss.deck_key != deck_key:
        build_deck(active_cats, shuffle, missed_only)
        ss.deck_key = deck_key


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
st.title("\U0001F9EA ChemCards")

ss = st.session_state

# A finish was requested last run; do the cookie write now, then the splash gate
# below ends this run with st.stop() (which flushes the pending cookie set to the
# browser -- unlike st.rerun(), which would discard it).
if ss.finish_pending:
    ss.finish_pending = False
    finish_daily()

# Completion splash takes over the page until dismissed.
if ss.show_splash is not None:
    render_splash(ss.show_splash)
    st.stop()

if not ss.deck:
    st.warning("No cards in the deck. Enable a group in the sidebar "
               "(or turn off 'Review missed only').")
    st.stop()

if ss.active_mode == "daily":
    _ans, _tot = len(ss.daily_answered), len(ss.deck)
    st.info(f"\U0001F525 Daily session \u2014 {_ans}/{_tot} answered")

card = COMPOUNDS[ss.deck[ss.pos]]
deck_idx = ss.deck[ss.pos]
struct_front = mode.startswith("Structure")
st.caption(f"Card {ss.pos + 1} / {len(ss.deck)}  \u2022  "
           f"{' \u00b7 '.join(card['categories'])}"
           f"  \u2022  attempts: {ss.attempts}")


# ===========================================================================
# MODE A: Structure -> Name  (type the answer)
# ===========================================================================
if struct_front:
    left, right = st.columns([1, 1])
    with left:
        st.image(render_structure(card["smiles"]))

    with right:
        if not ss.answered:
            with st.form(key=f"name_form_{ss.nonce}", clear_on_submit=False):
                guess_name = st.text_input("Name this compound",
                                           placeholder="e.g. dopamine, Phe, aspirin")
                checked = st.form_submit_button("Check", type="primary",
                                                use_container_width=True)
            if checked and guess_name.strip():
                ss.attempts += 1
                if _norm(guess_name) in card["aliases"]:
                    grade(True)
                    st.rerun()
                else:
                    ratio = difflib.SequenceMatcher(
                        None, _norm(guess_name), _norm(card["name"])).ratio()
                    if ratio > 0.8:
                        st.warning("So close \u2014 check your spelling and try again.")
                    else:
                        st.error("Not quite. Try again, or reveal the answer.")
            if st.button("Reveal answer / I don't know", use_container_width=True):
                grade(False)
                st.rerun()
        else:
            if ss.solved:
                st.success(f"Correct \u2014 **{card['name']}** "
                           f"(in {ss.attempts} attempt{'s' if ss.attempts != 1 else ''}).")
            else:
                st.error(f"The answer is **{card['name']}**.")
            if card.get("abbrev"):
                st.caption(card["abbrev"])
            if show_hints:
                st.markdown(f"**{card['formula']}**  \u2022  MW {card['mw']}")
            st.code(card["canonical"], language=None)

# ===========================================================================
# MODE B: Name -> Structure  (draw it in Ketcher)
# ===========================================================================
else:
    st.markdown(f"### {card['name']}")
    if card.get("abbrev"):
        st.caption(card["abbrev"])
    st.caption("Draw the structure, then click **Apply** in the editor.")
    # On narrow screens (mobile) the editor is wider than the viewport and its
    # atom palette gets clipped. Wrap it in a keyed container (stable CSS class
    # st-key-ketcher_wrap), pin a minimum width on the iframe, and let the
    # wrapper scroll horizontally. On desktop the natural width exceeds this, so
    # nothing changes there.
    st.markdown(
        """
        <style>
        .st-key-ketcher_wrap { overflow-x: auto !important;
                               -webkit-overflow-scrolling: touch; }
        .st-key-ketcher_wrap iframe { min-width: 720px !important; }
        </style>
        """,
        unsafe_allow_html=True,
    )
    with st.container(key="ketcher_wrap"):
        drawn = st_ketcher(key=f"ket_{ss.nonce}", height=520)
    # Process a drawing only once: Ketcher keeps returning the same SMILES on
    # every rerun, so guard on `last_processed` to avoid re-counting (and the
    # rerun loop that hid the feedback).
    if drawn and not ss.answered and drawn != ss.last_processed:
        gm = Chem.MolFromSmiles(drawn)
        if gm is None:
            st.warning("Couldn't parse that drawing \u2014 redraw and Apply again.")
        else:
            ss.last_processed = drawn
            ss.attempts += 1
            ss.last_guess = drawn
            if match_level(card["smiles"], gm) in ("exact", "constitution"):
                grade(True)
    if not ss.answered:
        if st.button("Reveal answer / I give up", use_container_width=True):
            grade(False)

    tm = Chem.MolFromSmiles(card["smiles"])

    # Feedback on the most recent drawn guess (closeness), if any.
    if ss.last_guess:
        gm = Chem.MolFromSmiles(ss.last_guess)
        lvl = match_level(card["smiles"], gm)
        if lvl == "exact":
            st.success("Exact match \u2014 structure and stereochemistry. \U0001F9EA")
        elif lvl == "constitution":
            st.success("Right constitution! (Stereochemistry differs from the "
                       "reference \u2014 fine for most flashcard purposes.)")
        else:
            st.info(f"Tanimoto to target (ECFP4): **{tanimoto(tm, gm):.2f}**")

        gf = rdMolDescriptors.CalcMolFormula(gm)
        st.markdown(f"Reference **{card['formula']}**  vs.  yours **{gf}**")
        st.dataframe(formula_table(tm, gm), hide_index=True,
                     use_container_width=True)
        if lvl == "no":
            for line in descriptor_deltas(tm, gm):
                if not line.startswith("\u2713"):
                    st.write("\u2022 " + line)
            st.caption("Atoms/bonds in the maximum common substructure are "
                       "highlighted on your drawing: :blue[blue] = tight match, "
                       ":red[red] = loose match (valence/aromaticity off).")
            st.image(render_comparison(card["smiles"], ss.last_guess))

    # Reveal the target once the card is answered.
    if ss.answered:
        st.divider()
        if ss.solved:
            st.success(f"Solved in {ss.attempts} attempt"
                       f"{'s' if ss.attempts != 1 else ''}.")
        else:
            st.error("Here's the reference structure:")
        st.image(render_structure(card["smiles"]))
        if show_hints:
            st.markdown(f"**{card['formula']}**  \u2022  MW {card['mw']}")
        st.code(card["canonical"], language=None)
    elif not ss.last_guess:
        st.caption("Draw the molecule above and click Apply to see how "
                   "close you are.")


# ---------------------------------------------------------------------------
# Navigation
# ---------------------------------------------------------------------------
st.divider()
n1, n2 = st.columns(2)
if n1.button("\u2190 Previous", use_container_width=True):
    advance(-1)
    st.rerun()
if n2.button("Next \u2192", use_container_width=True):
    advance(1)
    st.rerun()

if ss.active_mode == "daily":
    answered, total = len(ss.daily_answered), len(ss.deck)
    if answered >= total:
        if st.button("Finish today's session \U0001F389", type="primary",
                     use_container_width=True):
            # Don't write the cookie here: a set() followed by st.rerun() gets
            # discarded before it reaches the browser. Flag it and let the next
            # run do the write, then end on st.stop() (which flushes).
            ss.finish_pending = True
            st.rerun()
    else:
        st.caption(f"Answer all {total} cards to finish "
                   f"({total - answered} to go).")


def bind_enter_to_next(armed: bool):
    """Press Enter to click the 'Next' button, but only once a card is answered.

    The keydown listener is attached to the parent document a single time
    (guarded by a flag on the document); each rerun just re-sets whether it's
    'armed'. While unanswered, Enter is left alone so it still submits the
    name-entry form / does nothing in the drawing pane.
    """
    components.html(
        f"""
        <script>
        const doc = window.parent.document;
        if (!doc.__enterNextBound) {{
            doc.__enterNextBound = true;
            doc.addEventListener('keydown', function(e) {{
                if (e.key === 'Enter' && doc.__enterNextArmed) {{
                    for (const b of doc.querySelectorAll('button')) {{
                        if (b.innerText.trim().startsWith('Next')) {{
                            e.preventDefault();
                            b.click();
                            break;
                        }}
                    }}
                }}
            }});
        }}
        doc.__enterNextArmed = {str(armed).lower()};
        </script>
        """,
        height=0, width=0,
    )


bind_enter_to_next(ss.answered)
