"""Build an updated WHO seasonal vaccine-composition table.

Historical rows (through 2023-2024 NH / 2023 SH) are inherited from the
repository's curated CSV; 2024 SH onward are transcribed from WHO official
composition pages linked below. Output data stays under VAXSEER DATA.
"""
import argparse
import csv
import json
from pathlib import Path

WHO_HOME = "https://www.who.int/teams/global-influenza-programme/vaccines/who-recommendations"

NEW_ROWS = [
    {
        "Season": "2024 SH", "Hemisphere": "SH", "Recommendation_Date": "2023-09-29",
        "Platform": "Egg-based", "A/H1N1": "A/Victoria/4897/2022", "A/H3N2": "A/Thailand/8/2022",
        "B/Victoria": "B/Austria/1359417/2021", "B/Yamagata": "B/Phuket/3073/2013",
        "Valency_Note": "Trivalent recommended; B/Yamagata unchanged only where quadrivalent vaccines remain in use.",
        "Source_URL": "https://www.who.int/publications/m/item/recommended-composition-of-influenza-virus-vaccines-for-use-in-the-2024-southern-hemisphere-influenza-season",
        "Source_Basis": "WHO official report",
    },
    {
        "Season": "2024 SH", "Hemisphere": "SH", "Recommendation_Date": "2023-09-29",
        "Platform": "Cell/recombinant", "A/H1N1": "A/Wisconsin/67/2022", "A/H3N2": "A/Massachusetts/18/2022",
        "B/Victoria": "B/Austria/1359417/2021", "B/Yamagata": "B/Phuket/3073/2013",
        "Valency_Note": "Trivalent recommended; B/Yamagata unchanged only where quadrivalent vaccines remain in use.",
        "Source_URL": "https://www.who.int/publications/m/item/recommended-composition-of-influenza-virus-vaccines-for-use-in-the-2024-southern-hemisphere-influenza-season",
        "Source_Basis": "WHO official report",
    },
    {
        "Season": "2024-2025 NH", "Hemisphere": "NH", "Recommendation_Date": "2024-02-23",
        "Platform": "Egg-based", "A/H1N1": "A/Victoria/4897/2022", "A/H3N2": "A/Thailand/8/2022",
        "B/Victoria": "B/Austria/1359417/2021", "B/Yamagata": "B/Phuket/3073/2013",
        "Valency_Note": "Trivalent recommended; B/Yamagata unchanged only where quadrivalent vaccines remain in use.",
        "Source_URL": "https://www.who.int/publications/m/item/recommended-composition-of-influenza-virus-vaccines-for-use-in-the-2024-2025-northern-hemisphere-influenza-season",
        "Source_Basis": "WHO official report",
    },
    {
        "Season": "2024-2025 NH", "Hemisphere": "NH", "Recommendation_Date": "2024-02-23",
        "Platform": "Cell/recombinant", "A/H1N1": "A/Wisconsin/67/2022", "A/H3N2": "A/Massachusetts/18/2022",
        "B/Victoria": "B/Austria/1359417/2021", "B/Yamagata": "B/Phuket/3073/2013",
        "Valency_Note": "Trivalent recommended; B/Yamagata unchanged only where quadrivalent vaccines remain in use.",
        "Source_URL": "https://www.who.int/publications/m/item/recommended-composition-of-influenza-virus-vaccines-for-use-in-the-2024-2025-northern-hemisphere-influenza-season",
        "Source_Basis": "WHO official report",
    },
    {
        "Season": "2025 SH", "Hemisphere": "SH", "Recommendation_Date": "2024-09-27",
        "Platform": "Egg-based", "A/H1N1": "A/Victoria/4897/2022", "A/H3N2": "A/Croatia/10136RV/2023",
        "B/Victoria": "B/Austria/1359417/2021", "B/Yamagata": "B/Phuket/3073/2013",
        "Valency_Note": "Trivalent recommended; B/Yamagata unchanged only where quadrivalent vaccines remain in use.",
        "Source_URL": "https://www.who.int/publications/m/item/recommended-composition-of-influenza-virus-vaccines-for-use-in-the-2025-southern-hemisphere-influenza-season",
        "Source_Basis": "WHO official report",
    },
    {
        "Season": "2025 SH", "Hemisphere": "SH", "Recommendation_Date": "2024-09-27",
        "Platform": "Cell/recombinant/nucleic-acid", "A/H1N1": "A/Wisconsin/67/2022", "A/H3N2": "A/District of Columbia/27/2023",
        "B/Victoria": "B/Austria/1359417/2021", "B/Yamagata": "B/Phuket/3073/2013",
        "Valency_Note": "Trivalent recommended; B/Yamagata unchanged only where quadrivalent vaccines remain in use.",
        "Source_URL": "https://www.who.int/publications/m/item/recommended-composition-of-influenza-virus-vaccines-for-use-in-the-2025-southern-hemisphere-influenza-season",
        "Source_Basis": "WHO official report",
    },
    {
        "Season": "2025-2026 NH", "Hemisphere": "NH", "Recommendation_Date": "2025-02-28",
        "Platform": "Egg-based", "A/H1N1": "A/Victoria/4897/2022", "A/H3N2": "A/Croatia/10136RV/2023",
        "B/Victoria": "B/Austria/1359417/2021", "B/Yamagata": "B/Phuket/3073/2013",
        "Valency_Note": "Trivalent recommended; B/Yamagata unchanged only where quadrivalent vaccines remain in use.",
        "Source_URL": "https://www.who.int/publications/m/item/recommended-composition-of-influenza-virus-vaccines-for-use-in-the-2025-2026-nh-influenza-season",
        "Source_Basis": "WHO official report",
    },
    {
        "Season": "2025-2026 NH", "Hemisphere": "NH", "Recommendation_Date": "2025-02-28",
        "Platform": "Cell/recombinant/nucleic-acid", "A/H1N1": "A/Wisconsin/67/2022", "A/H3N2": "A/District of Columbia/27/2023",
        "B/Victoria": "B/Austria/1359417/2021", "B/Yamagata": "B/Phuket/3073/2013",
        "Valency_Note": "Trivalent recommended; B/Yamagata unchanged only where quadrivalent vaccines remain in use.",
        "Source_URL": "https://www.who.int/publications/m/item/recommended-composition-of-influenza-virus-vaccines-for-use-in-the-2025-2026-nh-influenza-season",
        "Source_Basis": "WHO official report",
    },
    {
        "Season": "2026 SH", "Hemisphere": "SH", "Recommendation_Date": "2025-09-26",
        "Platform": "Egg-based", "A/H1N1": "A/Missouri/11/2025", "A/H3N2": "A/Singapore/GP20238/2024",
        "B/Victoria": "B/Austria/1359417/2021", "B/Yamagata": "B/Phuket/3073/2013",
        "Valency_Note": "Trivalent recommended; B/Yamagata unchanged only where quadrivalent vaccines remain in use.",
        "Source_URL": "https://www.who.int/publications/m/item/recommended-composition-of-influenza-virus-vaccines-for-use-in-the-2026-southern-hemisphere-influenza-season",
        "Source_Basis": "WHO official report",
    },
    {
        "Season": "2026 SH", "Hemisphere": "SH", "Recommendation_Date": "2025-09-26",
        "Platform": "Cell/recombinant/nucleic-acid", "A/H1N1": "A/Missouri/11/2025", "A/H3N2": "A/Sydney/1359/2024",
        "B/Victoria": "B/Austria/1359417/2021", "B/Yamagata": "B/Phuket/3073/2013",
        "Valency_Note": "Trivalent recommended; B/Yamagata unchanged only where quadrivalent vaccines remain in use.",
        "Source_URL": "https://www.who.int/publications/m/item/recommended-composition-of-influenza-virus-vaccines-for-use-in-the-2026-southern-hemisphere-influenza-season",
        "Source_Basis": "WHO official report",
    },
    {
        "Season": "2026-2027 NH", "Hemisphere": "NH", "Recommendation_Date": "2026-02-27",
        "Platform": "Egg-based", "A/H1N1": "A/Missouri/11/2025", "A/H3N2": "A/Darwin/1454/2025",
        "B/Victoria": "B/Tokyo/EIS13-175/2025", "B/Yamagata": "B/Phuket/3073/2013",
        "Valency_Note": "Trivalent recommended; B/Yamagata unchanged only where quadrivalent vaccines remain in use.",
        "Source_URL": "https://www.who.int/publications/m/item/recommended-composition-of-influenza-virus-vaccines-for-use-in-the-2026-2027-northern-hemisphere-influenza-season",
        "Source_Basis": "WHO official report",
    },
    {
        "Season": "2026-2027 NH", "Hemisphere": "NH", "Recommendation_Date": "2026-02-27",
        "Platform": "Cell/recombinant/nucleic-acid", "A/H1N1": "A/Missouri/11/2025", "A/H3N2": "A/Darwin/1415/2025",
        "B/Victoria": "B/Pennsylvania/14/2025", "B/Yamagata": "B/Phuket/3073/2013",
        "Valency_Note": "Trivalent recommended; B/Yamagata unchanged only where quadrivalent vaccines remain in use.",
        "Source_URL": "https://www.who.int/publications/m/item/recommended-composition-of-influenza-virus-vaccines-for-use-in-the-2026-2027-northern-hemisphere-influenza-season",
        "Source_Basis": "WHO official report",
    },
    {
        "Season": "2027 SH", "Hemisphere": "SH", "Recommendation_Date": "2026-09-25",
        "Platform": "Egg-based", "A/H1N1": "A/Missouri/11/2025", "A/H3N2": "A/Darwin/1454/2025",
        "B/Victoria": "B/Tokyo/EIS13-175/2025", "B/Yamagata": "B/Phuket/3073/2013",
        "Valency_Note": "Trivalent recommended; B/Yamagata unchanged only where quadrivalent vaccines remain in use.",
        "Source_URL": "https://www.who.int/publications/m/item/recommended-composition-of-influenza-virus-vaccines-for-use-in-the-2027-southern-hemisphere-influenza-season",
        "Source_Basis": "WHO official report",
    },
    {
        "Season": "2027 SH", "Hemisphere": "SH", "Recommendation_Date": "2026-09-25",
        "Platform": "Cell/recombinant/nucleic-acid", "A/H1N1": "A/Missouri/11/2025", "A/H3N2": "A/Darwin/1415/2025",
        "B/Victoria": "B/Pennsylvania/14/2025", "B/Yamagata": "B/Phuket/3073/2013",
        "Valency_Note": "Trivalent recommended; B/Yamagata unchanged only where quadrivalent vaccines remain in use.",
        "Source_URL": "https://www.who.int/publications/m/item/recommended-composition-of-influenza-virus-vaccines-for-use-in-the-2027-southern-hemisphere-influenza-season",
        "Source_Basis": "WHO official report",
    },
]

FIELDS = ["Season", "Hemisphere", "Recommendation_Date", "Platform", "A/H1N1", "A/H3N2", "B/Victoria", "B/Yamagata", "Valency_Note", "Source_URL", "Source_Basis"]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--legacy-csv", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output = args.output_dir / "WHO_seasonal_vaccine_composition_2010_2027.csv"

    rows = []
    with args.legacy_csv.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.reader(f)
        headers = next(reader)
        columns = {name: index for index, name in enumerate(headers) if name}
        for raw in reader:
            season = (raw[columns["Season"]] or "").strip()
            h1 = (raw[columns["A/H1N1"]] or "").strip()
            if not season or not h1:
                continue
            hemi = "NH" if season.endswith("NH") else "SH"
            platform = (raw[5] if len(raw) > 5 else "").strip() or "Not specified in source table"
            rows.append({
                "Season": season, "Hemisphere": hemi, "Recommendation_Date": "",
                "Platform": platform, "A/H1N1": h1,
                "A/H3N2": (raw[columns["A/H3N2"]] or "").strip(),
                "B/Victoria": (raw[columns["B/Victoria"]] or "").strip(),
                "B/Yamagata": (raw[columns["B/Yamagata"]] or "").strip(),
                "Valency_Note": "Historical row retained from the repository's curated recommendation table; publication date and row-level WHO link not supplied there.",
                "Source_URL": WHO_HOME,
                "Source_Basis": "Repository curated historical table; official WHO season pages not individually reverified",
            })

    existing_seasons = {r["Season"] for r in rows}
    new_seasons = {r["Season"] for r in NEW_ROWS}
    duplicate_seasons = existing_seasons & new_seasons
    if duplicate_seasons:
        raise SystemExit("Supplement duplicates existing seasons: " + ", ".join(sorted(duplicate_seasons)))
    rows.extend(NEW_ROWS)
    with output.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)

    summary = {
        "as_of_date": "2026-09-28",
        "legacy_input": str(args.legacy_csv),
        "output_csv": str(output),
        "rows": len(rows),
        "legacy_rows": len(rows) - len(NEW_ROWS),
        "officially_added_rows": len(NEW_ROWS),
        "newest_northern_hemisphere_season": "2026-2027 NH",
        "newest_southern_hemisphere_season": "2027 SH",
        "historical_coverage_note": "The main multi-component table starts at 2010-2011 NH / 2010 SH. WHO archive H1N1 recommendations from 1998-2010 are supplied separately; archive rows include H1N1 only.",
        "archive_supplement_csv": str(args.output_dir / "WHO_H1N1_archive_1998_2010.csv"),
        "sequence_mapping_csv": str(args.output_dir / "WHO_H1N1_sequence_mapping_to_HIN1_MERGED.csv"),
        "official_archive_url": "https://www.who.int/teams/global-influenza-programme/vaccines/who-recommendations/recommendations-for-influenza-vaccine-composition-archive",
        "latest_recommendations_index": WHO_HOME,
        "b_yamagata_note": "Since WHO no longer considers a B/Yamagata antigen warranted for trivalent vaccines, B/Phuket/3073/2013 is retained only as the unchanged quadrivalent component where quadrivalent vaccines remain in use.",
    }
    (args.output_dir / "UPDATE_REPORT.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
