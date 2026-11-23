#!/usr/bin/env python3
"""GM12878 100 kb compartment knowledge graph.

A compact, indexed materialization of genomic-bin and compartment-membership nodes
and edges in SQLite. Hi-C INTERACTS_WITH edges remain in the hictile store and are
joined with compartment assignments on demand; avoids duplicating massive matrices.
Standard library only for building / label queries. hictile is required for `contacts`.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sqlite3
import sys
from pathlib import Path

RES = 100_000
NOT_LABELS = {"", "NA", "N/A", "NAN", "NONE", "NULL", ".", "-1"}


def fail(message):
    raise SystemExit("ERROR: " + message)


def clean_label(value):
    if value is None:
        return None
    value = str(value).strip()
    return value if value.upper() not in NOT_LABELS else None


def rows_tsv(path):
    with open(path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            fail(f"no tabular header in {path}")
        for row in reader:
            yield {k.strip(): v for k, v in row.items() if k is not None}


def node_id(row):
    key = "node" if "node" in row else "gid" if "gid" in row else None
    if key is None:
        fail("node column missing (expected 'node' or 'gid')")
    return int(row[key])


def initialize(conn):
    conn.executescript("""
        PRAGMA foreign_keys=ON;
        CREATE TABLE IF NOT EXISTS metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL);
        CREATE TABLE IF NOT EXISTS bins (
            node INTEGER PRIMARY KEY, chrom TEXT NOT NULL,
            start INTEGER NOT NULL, end INTEGER NOT NULL,
            UNIQUE(chrom, start), CHECK(start >= 0 AND end > start)
        );
        CREATE INDEX IF NOT EXISTS bins_by_region ON bins(chrom, start);
        CREATE TABLE IF NOT EXISTS compartments (
            method TEXT NOT NULL, label TEXT NOT NULL,
            PRIMARY KEY(method, label)
        );
        CREATE TABLE IF NOT EXISTS membership (
            node INTEGER NOT NULL REFERENCES bins(node),
            method TEXT NOT NULL, label TEXT NOT NULL,
            PRIMARY KEY(node, method),
            FOREIGN KEY(method, label) REFERENCES compartments(method, label)
        );
        CREATE INDEX IF NOT EXISTS membership_by_compartment ON membership(method, label, node);
    """)


def build(args):
    dest = Path(args.db)
    if dest.exists() and not args.overwrite:
        fail(f"{dest} already exists; use --overwrite to rebuild intentionally")
    files = sorted(Path(args.clusters_dir).glob("*_clusters.tsv"))
    if not files:
        fail(f"no *_clusters.tsv files in {args.clusters_dir}")
    if not Path(args.nodes).is_file():
        fail(f"nodes file not found: {args.nodes}")
    tmp = dest.with_name(dest.name + ".building")
    if tmp.exists():
        fail(f"temporary build exists: {tmp}; inspect/remove before retry")
    dest.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(tmp)
    counts = {"bins": 0, "memberships": 0, "files": len(files)}
    try:
        initialize(conn)
        conn.execute("INSERT INTO metadata VALUES (?,?)", ("cell_line", "GM12878"))
        conn.execute("INSERT INTO metadata VALUES (?,?)", ("resolution_bp", str(RES)))
        conn.execute("INSERT INTO metadata VALUES (?,?)", ("source_nodes", str(args.nodes)))
        known_nodes = {}
        for row in rows_tsv(args.nodes):
            node = node_id(row)
            chrom = row["chrom"].strip()
            start = int(row["start"])
            end = int(row["end"])
            if start % RES:
                fail(f"node {node} does not start on a 100kb boundary")
            if end <= start or end > start + RES:
                fail(f"node {node} has invalid coordinates: {start}..{end}")
            conn.execute("INSERT INTO bins VALUES (?,?,?,?)", (node, chrom, start, end))
            known_nodes[node] = (chrom, start, end)
            counts["bins"] += 1
        # Reference labels repeat in every clustering export: deduplicate and
        # fail on disagreements rather than silently overwriting them.
        assignments = {}
        for path in files:
            method = path.name.removesuffix("_clusters.tsv")
            for row in rows_tsv(path):
                node = node_id(row)
                if node not in known_nodes:
                    fail(f"node {node} from {path} missing from nodes index")
                if all(k in row for k in ("chrom", "start", "end")):
                    expected = known_nodes[node]
                    observed = (row["chrom"].strip(), int(row["start"]), int(row["end"]))
                    if expected != observed:
                        fail(f"node/coordinate disagreement for node {node}: {path}")
                for which, column in ((method, "cluster"), ("Rao", "rao"), ("SCI_published", "sci_pub")):
                    if column not in row:
                        if which == method:
                            fail(f"missing 'cluster' column: {path}")
                        continue
                    label = clean_label(row[column])
                    if label is None:
                        continue
                    key = (node, which)
                    previous = assignments.get(key)
                    if previous is not None and previous != label:
                        fail(f"conflicting compartment labels for node={node}, method={which}: {previous} vs {label}")
                    assignments[key] = label
        for method, label in sorted(set((m, lab) for (_n, m), lab in assignments.items())):
            conn.execute("INSERT INTO compartments VALUES (?,?)", (method, label))
        conn.executemany("INSERT INTO membership VALUES (?,?,?)",
                         ((node, method, label) for (node, method), label in assignments.items()))
        counts["memberships"] = len(assignments)
        conn.commit()
        conn.close()
        if dest.exists():
            dest.unlink()
        tmp.replace(dest)
        print(json.dumps({"status": "built", "database": str(dest), **counts}, indent=2))
    except BaseException:
        conn.close()
        if tmp.exists():
            tmp.unlink()
        raise


def connection(path):
    if not Path(path).is_file():
        fail(f"knowledge graph database not found: {path}")
    con = sqlite3.connect(f"file:{Path(path).resolve()}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    return con


def query_rows(con, sql, params=()):
    return [dict(row) for row in con.execute(sql, params)]


def emit(data):
    print(json.dumps(data, indent=2, ensure_ascii=False))


def parse_region(text, require_interval=False):
    m = re.fullmatch(r"(chr(?:[0-9]+|X|Y|M)):(\d+)(?:-(\d+))?", text, flags=re.I)
    if not m:
        fail("region must be chrN:POSITION or chrN:START-END (coordinates in bp)")
    chrom, x, y = m.group(1), int(m.group(2)), m.group(3)
    if require_interval and y is None:
        fail("expected chrN:START-END region")
    end = int(y) if y is not None else x + 1
    if end <= x:
        fail("region end must be greater than start")
    return chrom, x, end


def show_methods(args):
    with connection(args.db) as con:
        emit(query_rows(con, """SELECT c.method, COUNT(m.node) AS labeled_bins,
                COUNT(DISTINCT c.label) AS labels
                FROM compartments c LEFT JOIN membership m
                ON c.method=m.method AND c.label=m.label
                GROUP BY c.method ORDER BY c.method"""))


def summary(args):
    with connection(args.db) as con:
        emit(query_rows(con, """SELECT label, COUNT(*) AS bins
                 FROM membership WHERE method=?
                 GROUP BY label ORDER BY bins DESC, label""", (args.method,)))


def bin_query(args):
    chrom, start, end = parse_region(args.region)
    with connection(args.db) as con:
        bins = query_rows(con, """SELECT * FROM bins
                  WHERE chrom=? AND start < ? AND end > ? ORDER BY start LIMIT ?""",
                          (chrom, end, start, args.limit))
        for node in bins:
            node["compartments"] = {r["method"]: r["label"] for r in query_rows(
                con, "SELECT method, label FROM membership WHERE node=? ORDER BY method", (node["node"],))}
        emit(bins)


def compartment(args):
    sql = """SELECT b.node, b.chrom, b.start, b.end
             FROM membership AS m JOIN bins AS b ON b.node=m.node
             WHERE m.method=? AND m.label=?"""
    params = [args.method, args.label]
    if args.chrom:
        sql += " AND b.chrom=?"
        params.append(args.chrom)
    sql += " ORDER BY b.chrom, b.start LIMIT ?"
    params.append(args.limit)
    with connection(args.db) as con:
        count_sql = """SELECT COUNT(*) FROM membership m JOIN bins b ON b.node=m.node
                       WHERE m.method=? AND m.label=?"""
        count_params = [args.method, args.label]
        if args.chrom:
            count_sql += " AND b.chrom=?"
            count_params.append(args.chrom)
        emit({"method": args.method, "label": args.label,
              "count": con.execute(count_sql, count_params).fetchone()[0],
              "regions": query_rows(con, sql, params)})


def export_nodes_edges(args):
    """Export an interoperable property-graph schema without copying huge Hi-C edges."""
    folder = Path(args.out_dir)
    folder.mkdir(parents=True, exist_ok=True)
    definitions = [
        ("genomic_bins.csv", ["id", "chrom", "start", "end"],
         "SELECT 'bin:' || node AS id, chrom, start, end FROM bins ORDER BY node"),
        ("compartments.csv", ["id", "method", "label"],
         "SELECT 'compartment:' || method || ':' || label AS id, method, label FROM compartments ORDER BY method,label"),
        ("membership_edges.csv", ["source", "target", "relation"],
         "SELECT 'bin:' || node AS source, 'compartment:' || method || ':' || label AS target, "
         "'MEMBER_OF' AS relation FROM membership ORDER BY node, method"),
    ]
    stats = {}
    with connection(args.db) as con:
        for filename, cols, stmt in definitions:
            count = 0
            with (folder / filename).open('w', newline='', encoding='utf-8') as fh:
                w = csv.writer(fh)
                w.writerow(cols)
                for result in con.execute(stmt):
                    w.writerow(tuple(result))
                    count += 1
            stats[filename] = count
    emit({"destination": str(folder), "exported": stats,
          "note": "Hi-C contact relationships remain in hictile; no giant edge-list copied"})


def annotated_contacts(args):
    try:
        from hictile import HicStore, Engine
    except ImportError:
        fail("hictile isn't importable. Set PYTHONPATH=/path/to/unzipped/hictile")
    q = args.hql.strip()
    cmd = q.split(maxsplit=1)[0].upper() if q else ""
    if cmd not in {"CONTACTS", "TOPK", "NEIGHBORS"}:
        fail("contacts supports read-only CONTACTS, TOPK, NEIGHBORS queries")
    if re.search(r"\b(INTO|FORMAT|BOTH)\b", q, re.I):
        fail("output-writing syntax isn't permitted in compartment queries")
    m = re.search(r"\bRES\s+(\d+)(kb|mb|bp)?\b", q, re.I)
    if m:
        num = int(m.group(1)) * {"KB":1000, "MB":1000000, "BP":1, None:1}[m.group(2).upper() if m.group(2) else None]
        if num != RES:
            fail("compartment knowledge graph only supports 100 kb HQL queries")
    else:
        q += " RES 100kb"
    if cmd == "TOPK":
        k = re.match(r"TOPK\s+(\d+)\b", q, re.I)
        if not k or int(k.group(1)) > args.limit:
            fail("TOPK must specify K <= --limit")
    else:
        lm = re.search(r"\bLIMIT\s+(\d+)\b", q, re.I)
        if lm and int(lm.group(1)) > args.limit:
            fail("HQL LIMIT exceeds --limit")
        if not lm:
            q += f" LIMIT {args.limit}"
    store = HicStore(args.store, cache_mb=256, threads=2)
    result = Engine(store).run(q)
    if result.kind != "table" or result.table is None:
        fail("HQL query must return a table")
    table = result.table.head(args.limit)
    if "bin1" not in table.columns or "bin2" not in table.columns:
        fail("HQL returned no bin1/bin2 columns")
    methods = args.method or ["Rao", "SCI_published"]
    nodes = set(int(x) for col in ("bin1", "bin2") for x in table[col])
    labels = {}
    node_meta = {}
    with connection(args.db) as con:
        node_list = sorted(nodes)
        for pos in range(0, len(node_list), 800):
            part = node_list[pos:pos+800]
            qmarks = ",".join("?" for _ in part)
            for r in con.execute(f"SELECT node, method, label FROM membership WHERE node IN ({qmarks})", part):
                labels[(int(r["node"]), r["method"])] = r["label"]
            for r in con.execute(f"SELECT node, chrom, start FROM bins WHERE node IN ({qmarks})", part):
                node_meta[int(r["node"])] = (r["chrom"], int(r["start"]))
    # A different chromosome order or bin-numbering convention in hictile would
    # silently attach biologically WRONG compartment labels. Reject such stores.
    for row in table.to_dict("records"):
        for idx in (1, 2):
            node = int(row[f"bin{idx}"])
            expected = node_meta.get(node)
            observed = (str(row[f"chrom{idx}"]), int(row[f"start{idx}"]))
            if expected != observed:
                fail(f"hictile store / KG node-ID mismatch: bin {node}, "
                     f"store={observed}, KG={expected}. Check chromosome order and bin rule.")
    data = []
    for row in table.to_dict("records"):
        a,b = int(row["bin1"]), int(row["bin2"])
        entry = {k: (v.item() if hasattr(v, "item") else v) for k,v in row.items()}
        for method in methods:
            key = re.sub(r"[^a-zA-Z0-9]+", "_", method)
            entry[f"{key}_compartment_1"] = labels.get((a, method))
            entry[f"{key}_compartment_2"] = labels.get((b, method))
            entry[f"{key}_same_compartment"] = (
                entry[f"{key}_compartment_1"] == entry[f"{key}_compartment_2"]
                if entry[f"{key}_compartment_1"] is not None and entry[f"{key}_compartment_2"] is not None else None
            )
        data.append(entry)
    if args.filter_method or args.label1 or args.label2 or args.only_same or args.only_different:
        scheme = args.filter_method or (methods[0] if len(methods) == 1 else None)
        if not scheme:
            fail("specify --filter-method to select the compartment-label scheme")
        if scheme not in methods:
            fail("--filter-method must also appear as one of the --method values")
        key = re.sub(r"[^a-zA-Z0-9]+", "_", scheme)
        selected = []
        for entry in data:
            a = entry[f"{key}_compartment_1"]
            b = entry[f"{key}_compartment_2"]
            if args.label1 and a != args.label1:
                continue
            if args.label2 and b != args.label2:
                continue
            if args.only_same and not (a is not None and b is not None and a == b):
                continue
            if args.only_different and not (a is not None and b is not None and a != b):
                continue
            selected.append(entry)
        data = selected
    emit({"cell_line": "GM12878", "resolution_bp": RES,
          "hql": q, "returned": len(data),
          "note": "Compartment filters apply to the HQL result set, NOT to all contacts in the full genome.",
          "contacts": data})


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("build", help="Build indexed compartment KG from node and cluster tables")
    p.add_argument("--nodes", required=True)
    p.add_argument("--clusters-dir", required=True)
    p.add_argument("--db", required=True)
    p.add_argument("--overwrite", action="store_true")
    p.set_defaults(func=build)
    for name,func in (("methods",show_methods),("summary",summary),("bin",bin_query),
                      ("compartment",compartment),("contacts",annotated_contacts),
                      ("export",export_nodes_edges)):
        p = sub.add_parser(name)
        p.add_argument("--db", required=True)
        if name in {"summary", "compartment"}:
            p.add_argument("--method", required=True)
        if name == "bin":
            p.add_argument("--region", required=True, help="chr1:100000 or chr1:100000-500000")
            p.add_argument("--limit", type=int, default=20)
        if name == "compartment":
            p.add_argument("--label", required=True)
            p.add_argument("--chrom")
            p.add_argument("--limit", type=int, default=20)
        if name == "export":
            p.add_argument("--out-dir", required=True)
        if name == "contacts":
            p.add_argument("--store", required=True, help="hictile .htl path")
            p.add_argument("--hql", required=True)
            p.add_argument("--method", action="append", help="membership scheme; repeatable")
            p.add_argument("--limit", type=int, default=100)
            p.add_argument("--filter-method", help="method to use for compartment filters")
            p.add_argument("--label1", help="only contacts whose first bin has this label")
            p.add_argument("--label2", help="only contacts whose second bin has this label")
            filters = p.add_mutually_exclusive_group()
            filters.add_argument("--only-same", action="store_true", help="within-compartment contacts")
            filters.add_argument("--only-different", action="store_true", help="between-compartment contacts")
        p.set_defaults(func=func)
    args = parser.parse_args(argv)
    if hasattr(args, "limit") and (args.limit < 1 or args.limit > 1000):
        fail("--limit must be 1 to 1000")
    args.func(args)


if __name__ == "__main__":
    main()
