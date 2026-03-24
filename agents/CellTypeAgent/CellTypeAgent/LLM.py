""" This file contains the code for calling all LLM APIs. Ref: https://github.com/snap-stanford/BioDiscoveryAgent/blob/master/LLM.py"""
import sys
import os
import time
import tiktoken
from run_metrics import get_active_run_tracker

print("Current working directory:", os.getcwd())


def _get_encoder(model_name):
    try:
        return tiktoken.encoding_for_model(model_name)
    except Exception:
        return tiktoken.get_encoding("cl100k_base")

try:
    import anthropic
    _anthropic_human_prompt = getattr(anthropic, "HUMAN_PROMPT", "\n\nHuman:")
    # Please setup your Anthropic API key
    claude_key_path = "APIs/claude_api_key.txt"
    print(f"Looking for Claude API key at: {os.path.abspath(claude_key_path)}")
    anthropic_client = anthropic.Client(api_key=open(claude_key_path).read().strip())
except Exception as e:
    anthropic = None
    _anthropic_human_prompt = "\n\nHuman:"
    anthropic_client = None
    print(e)
    print("Could not load anthropic API key.")

try:
    import openai
    from openai import OpenAI
    import yaml
    # Please setup your OpenAI API key
    try:
        openai_key_path = "APIs/openai_api_key.txt"
        print(f"Looking for OpenAI API key at: {os.path.abspath(openai_key_path)}")
        openai_client = OpenAI(
            api_key = open(openai_key_path).read().strip(),
            timeout=60*10, # o1-preview need a long time to respond
        )
    except FileNotFoundError:
        print(f"Error: OpenAI API key file not found at {os.path.abspath(openai_key_path)}")
        print("Please create the file and add your OpenAI API key")
        # Also check for the alternative filename mentioned in the error message
        alt_key_path = "APIs/openai_key.txt"
        print(f"Also checking alternative path: {os.path.abspath(alt_key_path)}")
        if os.path.exists(alt_key_path):
            print(f"Found key at alternative path: {alt_key_path}")
            print("There is a mismatch between the error message and the actual file being searched.")
        sys.exit(1)
except Exception as e:
    openai_client = None
    print(e)
    print("Could not initialize OpenAI client.")

try:
    import openai
    import yaml
    # Please setup your DeepSeek API key
    deepseek_key_path = "APIs/deepseek_api_key.yaml"
    print(f"Looking for DeepSeek API key at: {os.path.abspath(deepseek_key_path)}")
    key = yaml.safe_load(open(deepseek_key_path).read())
    deepseek_client = OpenAI(
        base_url = key["base_url"],
        api_key = key["api_key"]
    )
except Exception as e:
    deepseek_client = None
    print(e)
    print(f"Could not load deepseek API key from {deepseek_key_path}")


class TooLongPromptError(Exception):
    """Exception raised for errors in the prompt length."""
    def __init__(self, message="The prompt is too long."):
        self.message = message
        super().__init__(self.message)

class LLMError(Exception):
    """Exception raised for errors related to the LLM."""
    def __init__(self, message="An error occurred with the LLM."):
        self.message = message
        super().__init__(self.message)


def _safe_int(value, default=0):
    try:
        return int(value)
    except Exception:
        return default


def _get_usage_value(usage, key):
    if usage is None:
        return 0
    if isinstance(usage, dict):
        return _safe_int(usage.get(key, 0), 0)
    return _safe_int(getattr(usage, key, 0), 0)


def _usage_from_anthropic(message):
    usage = getattr(message, "usage", None)
    input_tokens = _get_usage_value(usage, "input_tokens")
    output_tokens = _get_usage_value(usage, "output_tokens")
    return input_tokens, output_tokens


def _usage_from_openai(response):
    usage = getattr(response, "usage", None)
    input_tokens = _get_usage_value(usage, "prompt_tokens")
    output_tokens = _get_usage_value(usage, "completion_tokens")
    return input_tokens, output_tokens


def _estimate_tokens(prompt, completion, model):
    encoder = _get_encoder(model)
    num_prompt_tokens = len(encoder.encode(prompt or ""))
    num_sample_tokens = len(encoder.encode(completion or ""))
    return num_prompt_tokens, num_sample_tokens


def _resolve_tokens(prompt, completion, model, input_tokens, output_tokens):
    if (input_tokens or 0) > 0 or (output_tokens or 0) > 0:
        return input_tokens, output_tokens, "api"
    est_input, est_output = _estimate_tokens(prompt, completion, model)
    return est_input, est_output, "estimated"


def _record_llm_metrics(provider, model, input_tokens, output_tokens, latency_seconds, log_file=None):
    tracker = get_active_run_tracker()
    if tracker is not None:
        tracker.record_llm_call(
            provider=provider,
            model=model,
            input_tokens=input_tokens,
            output_tokens=output_tokens,
            latency_seconds=latency_seconds,
            log_file=log_file,
            call_started_unix=time.time() - latency_seconds,
        )


def log_to_file(log_file, prompt, completion, model, input_tokens=None, output_tokens=None, latency_seconds=None, token_source="api"):
    """ Log the prompt and completion to a file."""
    if input_tokens is None or output_tokens is None:
        input_tokens, output_tokens = _estimate_tokens(prompt, completion, model)
        token_source = "estimated"

    with open(log_file, "a") as f:
        f.write("\n===================prompt=====================\n")
        f.write(prompt)
        f.write(f"\n==================={model} response =====================\n")
        f.write(completion)
        f.write("\n===================tokens=====================\n")
        f.write(f"Token source: {token_source}\n")
        f.write(f"Number of prompt tokens: {input_tokens}\n")
        f.write(f"Number of sampled tokens: {output_tokens}\n")
        f.write(f"Number of total tokens: {input_tokens + output_tokens}\n")
        if latency_seconds is not None:
            f.write(f"Latency seconds: {latency_seconds:.4f}\n")
        f.write("\n\n")


def complete_text_claude(prompt, stop_sequences=None, model="claude-3.5-sonnet", max_tokens_to_sample = 2000, temperature=0.5, log_file=None, **kwargs):
    """ Call the Claude API to complete a prompt."""

    if anthropic_client is None:
        raise LLMError("Anthropic client is not initialized.")

    if stop_sequences is None:
        stop_sequences = [_anthropic_human_prompt]

    start_time = time.perf_counter()
    try:
        message = anthropic_client.messages.create(
            model=model,
            max_tokens=max_tokens_to_sample,
            temperature=temperature,
            messages=[{"role": "user", "content": prompt}],
            stop_sequences=stop_sequences,
            **kwargs
        )
    except anthropic.APIStatusError as e:
        print(e)
        sys.exit()
        raise TooLongPromptError()
    except Exception as e:
        raise LLMError(str(e))

    latency_seconds = time.perf_counter() - start_time
    completion = message.content[0].text
    input_tokens, output_tokens = _usage_from_anthropic(message)
    input_tokens, output_tokens, token_source = _resolve_tokens(prompt, completion, model, input_tokens, output_tokens)

    if log_file is not None:
        log_to_file(
            log_file,
            prompt,
            completion,
            model,
            input_tokens=input_tokens,
            output_tokens=output_tokens,
            latency_seconds=latency_seconds,
            token_source=token_source,
        )
    _record_llm_metrics("anthropic", model, input_tokens, output_tokens, latency_seconds, log_file)
    return completion

def complete_text_openai(prompt, model="gpt-4o-mini-2024-07-18", temperature=0.5, log_file=None, **kwargs):

    """ Call the OpenAI API to complete a prompt."""
    if openai_client is None:
        raise LLMError("OpenAI client is not initialized.")

    raw_request = {
        "model": model,
        # "max_tokens": max_tokens_to_sample,
        **kwargs
    }
    
    # Add temperature only if supported by the model
    if model.startswith("o1"):
        temperature = 1.0 # o1 models only support temperature 1.0
        raw_request["temperature"] = temperature
    elif not model.startswith("o3-mini"):  # o3-mini models don't support temperature
        raw_request["temperature"] = temperature

    start_time = time.perf_counter()
    if model.startswith("gpt-3.5") or model.startswith("gpt-4") or model.startswith("gpt-5") or model.startswith("o1") or model.startswith("o3"):
        messages = [{"role": "user", "content": prompt}]
        response = openai_client.chat.completions.create(**{"messages": messages,**raw_request})
        completion = response.choices[0].message.content
    else:
        response = openai_client.chat.completions.create(model, **{"prompt": prompt,**raw_request})
        completion = response.choices[0].text

    latency_seconds = time.perf_counter() - start_time
    input_tokens, output_tokens = _usage_from_openai(response)
    input_tokens, output_tokens, token_source = _resolve_tokens(prompt, completion, model, input_tokens, output_tokens)

    if log_file is not None:
        log_to_file(
            log_file,
            prompt,
            completion,
            model,
            input_tokens=input_tokens,
            output_tokens=output_tokens,
            latency_seconds=latency_seconds,
            token_source=token_source,
        )
    _record_llm_metrics("openai", model, input_tokens, output_tokens, latency_seconds, log_file)
    return completion

def complete_text_deepseek(prompt, model="deepseek-r1", temperature=0.5, log_file=None, stream=False, **kwargs):
    """ Call the DeepSeek API to complete a prompt."""
    if deepseek_client is None:
        raise LLMError("DeepSeek client is not initialized.")

    max_attempts = 10
    for attempt in range(max_attempts):
        start_time = time.perf_counter()
        response = deepseek_client.chat.completions.create(model=f"deepseek-ai/{model}", messages=[{"role": "user", "content": prompt}], temperature=temperature, stream=stream, **kwargs)
        try:
            completion = response.choices[0].message.content
            latency_seconds = time.perf_counter() - start_time
            input_tokens, output_tokens = _usage_from_openai(response)
            input_tokens, output_tokens, token_source = _resolve_tokens(prompt, completion, model, input_tokens, output_tokens)

            if log_file is not None:
                log_to_file(
                    log_file,
                    prompt,
                    completion,
                    model,
                    input_tokens=input_tokens,
                    output_tokens=output_tokens,
                    latency_seconds=latency_seconds,
                    token_source=token_source,
                )
            _record_llm_metrics("deepseek", model, input_tokens, output_tokens, latency_seconds, log_file)
            return completion
        except Exception as e:
            print(e)
            print(f"Attempt {attempt+1} failed. Retrying...")
    raise LLMError("Failed to complete the prompt.")

def complete_text(prompt, log_file, model, **kwargs):
    """ Complete text using the specified model with appropriate API. """

    if model.startswith("claude"):
        completion = complete_text_claude(prompt, stop_sequences=[_anthropic_human_prompt, "Observation:"], log_file=log_file, model=model, **kwargs)
    elif model.startswith("gpt") or model.startswith("o1") or model.startswith("o3"):
        completion = complete_text_openai(prompt, log_file=log_file, model=model, **kwargs)
    elif model.startswith("deepseek"):
        completion = complete_text_deepseek(prompt, log_file=log_file, model=model, **kwargs)
        print(completion)
    return completion
