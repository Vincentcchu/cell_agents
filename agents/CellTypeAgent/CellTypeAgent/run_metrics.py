import csv
import json
import os
import time
from dataclasses import dataclass, asdict
from typing import Any, Dict, List, Optional


INPUT_COST_PER_MILLION = 1.1
OUTPUT_COST_PER_MILLION = 4.4


@dataclass
class LLMCallRecord:
    provider: str
    model: str
    input_tokens: int
    output_tokens: int
    total_tokens: int
    latency_seconds: float
    log_file: Optional[str] = None
    call_started_unix: Optional[float] = None


class RunMetricsTracker:
    def __init__(
        self,
        run_name: str,
        model: str,
        metadata: Optional[Dict[str, Any]] = None,
        input_cost_per_million: float = INPUT_COST_PER_MILLION,
        output_cost_per_million: float = OUTPUT_COST_PER_MILLION,
    ):
        self.run_name = run_name
        self.model = model
        self.metadata = metadata or {}
        self.input_cost_per_million = input_cost_per_million
        self.output_cost_per_million = output_cost_per_million
        self.started_at_unix = time.time()
        self._started_at_perf = time.perf_counter()
        self.llm_calls: List[LLMCallRecord] = []

    def record_llm_call(
        self,
        provider: str,
        model: str,
        input_tokens: int,
        output_tokens: int,
        latency_seconds: float,
        log_file: Optional[str] = None,
        call_started_unix: Optional[float] = None,
    ) -> None:
        input_tokens = int(input_tokens or 0)
        output_tokens = int(output_tokens or 0)
        self.llm_calls.append(
            LLMCallRecord(
                provider=provider,
                model=model,
                input_tokens=input_tokens,
                output_tokens=output_tokens,
                total_tokens=input_tokens + output_tokens,
                latency_seconds=float(latency_seconds or 0.0),
                log_file=log_file,
                call_started_unix=call_started_unix,
            )
        )

    def _aggregate(self) -> Dict[str, Any]:
        input_tokens = sum(call.input_tokens for call in self.llm_calls)
        output_tokens = sum(call.output_tokens for call in self.llm_calls)
        total_tokens = input_tokens + output_tokens
        annotation_time_seconds = sum(call.latency_seconds for call in self.llm_calls)
        total_time_seconds = time.perf_counter() - self._started_at_perf
        preprocessing_time_seconds = max(total_time_seconds - annotation_time_seconds, 0.0)
        input_cost_usd = (input_tokens / 1_000_000) * self.input_cost_per_million
        output_cost_usd = (output_tokens / 1_000_000) * self.output_cost_per_million
        total_cost_usd = input_cost_usd + output_cost_usd

        return {
            "run_name": self.run_name,
            "model": self.model,
            "started_at_unix": self.started_at_unix,
            "started_at_iso": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(self.started_at_unix)),
            "num_llm_calls": len(self.llm_calls),
            "input_tokens": input_tokens,
            "output_tokens": output_tokens,
            "total_tokens": total_tokens,
            "annotation_time_seconds": round(annotation_time_seconds, 6),
            "preprocessing_time_seconds": round(preprocessing_time_seconds, 6),
            "total_time_seconds": round(total_time_seconds, 6),
            "input_cost_per_million": self.input_cost_per_million,
            "output_cost_per_million": self.output_cost_per_million,
            "input_cost_usd": round(input_cost_usd, 8),
            "output_cost_usd": round(output_cost_usd, 8),
            "total_cost_usd": round(total_cost_usd, 8),
            "metadata": self.metadata,
            "llm_calls": [asdict(call) for call in self.llm_calls],
        }

    def finalize(self) -> Dict[str, Any]:
        return self._aggregate()

    def save_json(self, file_path: str) -> Dict[str, Any]:
        result = self.finalize()
        os.makedirs(os.path.dirname(file_path), exist_ok=True)
        with open(file_path, "w", encoding="utf-8") as f:
            json.dump(result, f, indent=2)
        return result

    def append_summary_csv(self, file_path: str, finalized: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        summary = finalized or self.finalize()
        flat = {
            "run_name": summary["run_name"],
            "model": summary["model"],
            "started_at_iso": summary["started_at_iso"],
            "num_llm_calls": summary["num_llm_calls"],
            "input_tokens": summary["input_tokens"],
            "output_tokens": summary["output_tokens"],
            "total_tokens": summary["total_tokens"],
            "annotation_time_seconds": summary["annotation_time_seconds"],
            "preprocessing_time_seconds": summary["preprocessing_time_seconds"],
            "total_time_seconds": summary["total_time_seconds"],
            "input_cost_usd": summary["input_cost_usd"],
            "output_cost_usd": summary["output_cost_usd"],
            "total_cost_usd": summary["total_cost_usd"],
        }

        os.makedirs(os.path.dirname(file_path), exist_ok=True)
        file_exists = os.path.exists(file_path)
        with open(file_path, "a", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(flat.keys()))
            if not file_exists:
                writer.writeheader()
            writer.writerow(flat)
        return flat


_ACTIVE_RUN_TRACKER: Optional[RunMetricsTracker] = None


def set_active_run_tracker(tracker: Optional[RunMetricsTracker]) -> None:
    global _ACTIVE_RUN_TRACKER
    _ACTIVE_RUN_TRACKER = tracker


def get_active_run_tracker() -> Optional[RunMetricsTracker]:
    return _ACTIVE_RUN_TRACKER


def clear_active_run_tracker() -> None:
    set_active_run_tracker(None)
