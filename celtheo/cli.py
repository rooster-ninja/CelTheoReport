"""Command line: python -m celtheo run [options]"""
import argparse
import sys


def main(argv=None):
    p = argparse.ArgumentParser(prog="celtheo", description="CelTheo occultation pipeline")
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="compute all events and write out/results.json + out/results.csv")
    r.add_argument("--csv", default="data/occultations.csv")
    r.add_argument("--out", default="out")
    r.add_argument("--offline", action="store_true", help="no network: use cached stars/soundings only")
    r.add_argument("--only", nargs="*", help="substring(s) of event ids to run, e.g. pikes-peak")
    r.add_argument("--no-tests", action="store_true", help="skip the ϵ/♈︎ test rows")
    r.add_argument("--refresh-stars", action="store_true", help="re-query SIMBAD for cached stars that lack proper motion")
    r.add_argument("--no-profiles", action="store_true", help="omit ray path / atmosphere profiles from JSON")
    ls = sub.add_parser("list", help="list event ids and warnings from the CSV")
    ls.add_argument("--csv", default="data/occultations.csv")
    a = p.parse_args(argv)

    if a.cmd == "list":
        from .data import load_events
        for e in load_events(a.csv):
            print(f"{e.id:60s} {e.utc:%Y-%m-%d %H:%M:%SZ}" if e.utc else e.id, "|", "; ".join(e.warnings))
        return
    from .api import run
    doc = run(a.csv, a.out, offline=a.offline, only=a.only, include_tests=not a.no_tests, profiles=not a.no_profiles, refresh_stars=a.refresh_stars)
    n = len(doc["events"])
    bad = sum(1 for e in doc["events"] if e.get("flags", {}).get("warnings"))
    print(f"\n{n} events -> {a.out}/results.json, {a.out}/results.csv ({bad} with warnings, {doc['runtime_s']} s)")


if __name__ == "__main__":
    sys.exit(main())
