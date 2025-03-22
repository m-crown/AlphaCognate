import { useState } from "react";
import { Combobox, Loader, TextInput, useCombobox } from '@mantine/core';
import { useDebouncedValue } from "@mantine/hooks";
import { useEffect } from "react";

async function fetchSearchResults(query: string) {
  const response = await fetch(`http://localhost:8000/search/?query=${encodeURIComponent(query)}`);
  if (!response.ok) throw new Error("Network error");
  return response.json();
}

export function AutoComplete() {
  const combobox = useCombobox();
  const [value, setValue] = useState("");
  const [data, setData] = useState<string[]>([]);
  const [loading, setLoading] = useState(false);
  const [empty, setEmpty] = useState(false);
  
  // Debounce the search value by 300ms
  const [debouncedValue] = useDebouncedValue(value, 300);

  // Fetch search results when debounced value changes
  useEffect(() => {
    if (debouncedValue.length > 1) {
      setLoading(true);
      fetchSearchResults(debouncedValue)
        .then((result) => {
          setData(result);
          setLoading(false);
          setEmpty(result.length === 0);
        })
        .catch(() => setLoading(false));
    }
  }, [debouncedValue]);

  return (
    <Combobox store={combobox} onOptionSubmit={(optionValue) => { setValue(optionValue); combobox.closeDropdown(); }}>
      <Combobox.Target>
        <TextInput
          label="Search Structures"
          placeholder="Type to search..."
          value={value}
          onChange={(event) => setValue(event.currentTarget.value)}
          onFocus={() => combobox.openDropdown()}
          rightSection={loading && <Loader size={18} />}
        />
      </Combobox.Target>

      <Combobox.Dropdown hidden={data.length === 0}>
        <Combobox.Options>
          {data.map((item) => <Combobox.Option value={item} key={item}>{item}</Combobox.Option>)}
          {empty && <Combobox.Empty>No results found</Combobox.Empty>}
        </Combobox.Options>
      </Combobox.Dropdown>
    </Combobox>
  );
}

export default AutoComplete;